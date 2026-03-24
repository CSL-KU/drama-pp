#!/usr/bin/env python3
"""
GF(2) DRAM bank-mapping solver from discovered same-bank address sets.

The legacy DRAMA++ workflow assumed one global XOR mapping for the entire
address space. This version keeps that path intact for symmetric systems and
adds region-aware solving for asymmetric DIMM layouts.
"""

import argparse
import itertools
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass
class BankSet:
    name: str
    addrs: List[int]
    rows: List[List[int]]


@dataclass(frozen=True)
class SelectorTerm:
    expected: int
    bits: Tuple[int, ...]

    def parity(self, addr: int) -> int:
        return sum((addr >> bit) & 1 for bit in self.bits) & 1

    def matches(self, addr: int) -> bool:
        return self.parity(addr) == self.expected

    def format_line(self) -> str:
        return 'selector {} {}'.format(self.expected, ' '.join(str(bit) for bit in self.bits))

    def describe(self) -> str:
        joined = ' ^ '.join(f'a{bit}' for bit in self.bits)
        return f'{joined}={self.expected}'


@dataclass(frozen=True)
class RegionSelector:
    terms: Tuple[SelectorTerm, ...]

    def matches(self, addr: int) -> bool:
        return all(term.matches(addr) for term in self.terms)

    def format_block(self) -> List[str]:
        lines = ['region']
        for term in self.terms:
            lines.append(term.format_line())
        return lines

    def describe(self) -> str:
        if not self.terms:
            return 'global'
        return ', '.join(term.describe() for term in self.terms)

    def complexity(self) -> int:
        return sum(len(term.bits) for term in self.terms)

    @staticmethod
    def from_mask_value(mask: int, value: int) -> 'RegionSelector':
        terms: List[SelectorTerm] = []
        bit = 0
        current_mask = mask
        while current_mask:
            if current_mask & 1:
                terms.append(SelectorTerm(expected=(value >> bit) & 1, bits=(bit,)))
            current_mask >>= 1
            bit += 1
        return RegionSelector(tuple(terms))


@dataclass
class RegionSolveResult:
    selector: Optional[RegionSelector]
    bank_sets: List[BankSet]
    matrix: List[List[int]]
    bank_codes: List[Tuple[int, ...]]
    mismatches: int
    total_addresses: int

    @property
    def bit_count(self) -> int:
        return len(self.matrix[0]) if self.matrix and self.matrix[0] else 0


@dataclass
class RegionProgress:
    result: RegionSolveResult
    unique_bank_codes: int
    unresolved_sets: int


@dataclass
class DiagnosticSplit:
    bits: Tuple[int, ...]
    child_sizes: Tuple[int, int]
    child_uniques: Tuple[int, int]


@dataclass
class BranchDiagnostic:
    selector: RegionSelector
    bank_set_names: List[str]
    total_sets: int
    unique_bank_codes: int
    reason: str
    suggested_splits: List[DiagnosticSplit]


@dataclass
class AutoRegionCandidate:
    results: List[RegionSolveResult]
    score: Tuple[int, int, int, int]
    diagnostics: List[BranchDiagnostic]


def parse_addr(text: str) -> Optional[int]:
    text = text.strip()
    if not text or text.startswith('#'):
        return None
    try:
        if text.startswith(('0x', '0X')):
            return int(text, 16)
        if any(ch in text.lower() for ch in 'abcdef'):
            return int(text, 16)
        return int(text, 10)
    except ValueError:
        return None


def read_addr_file(path: str) -> List[int]:
    addrs: List[int] = []
    with open(path, 'r', encoding='utf-8') as handle:
        for line in handle:
            parsed = parse_addr(line)
            if parsed is not None:
                addrs.append(parsed)
    if not addrs:
        raise ValueError(f'No valid addresses in {path}')
    return addrs


def build_bitvec(addr: int, lowbit: int, highbit: int) -> List[int]:
    return [(addr >> bit) & 1 for bit in range(lowbit, highbit + 1)]


def numeric_sort_key(path_text: str) -> Tuple[int, str]:
    stem = Path(path_text).stem
    digits = ''.join(ch for ch in stem if ch.isdigit())
    return (int(digits) if digits else sys.maxsize, path_text)


def load_bank_sets(file_paths: Sequence[str], lowbit: int, highbit: int) -> List[BankSet]:
    bank_sets: List[BankSet] = []
    for path in sorted(file_paths, key=numeric_sort_key):
        addrs = read_addr_file(path)
        rows = [build_bitvec(addr, lowbit, highbit) for addr in addrs]
        bank_sets.append(BankSet(name=Path(path).name, addrs=addrs, rows=rows))
    return bank_sets


def mat_rank_gf2(matrix: List[List[int]]) -> int:
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    n_rows = len(work)
    n_cols = len(work[0])
    rank = 0
    for col in range(n_cols):
        pivot = None
        for row in range(rank, n_rows):
            if work[row][col]:
                pivot = row
                break
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for row in range(n_rows):
            if row != rank and work[row][col]:
                for inner in range(col, n_cols):
                    work[row][inner] ^= work[rank][inner]
        rank += 1
        if rank == n_rows:
            break
    return rank


def nullspace_basis_gf2(diff_matrix: List[List[int]], nbits: int) -> List[List[int]]:
    if not diff_matrix:
        return [[1 if row == col else 0 for row in range(nbits)] for col in range(nbits)]

    work = [row[:] for row in diff_matrix]
    n_rows = len(work)
    pivot_cols: List[int] = []
    row = 0
    for col in range(nbits):
        pivot = None
        for candidate in range(row, n_rows):
            if work[candidate][col]:
                pivot = candidate
                break
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        pivot_cols.append(col)
        for candidate in range(n_rows):
            if candidate != row and work[candidate][col]:
                for inner in range(col, nbits):
                    work[candidate][inner] ^= work[row][inner]
        row += 1
        if row == n_rows:
            break

    pivot_set = set(pivot_cols)
    free_cols = [col for col in range(nbits) if col not in pivot_set]
    basis: List[List[int]] = []
    for free_col in free_cols:
        vector = [0] * nbits
        vector[free_col] = 1
        for pivot_row in range(len(pivot_cols) - 1, -1, -1):
            pivot_col = pivot_cols[pivot_row]
            accum = 0
            for col in range(pivot_col + 1, nbits):
                if work[pivot_row][col] and vector[col]:
                    accum ^= 1
            vector[pivot_col] = accum
        basis.append(vector)
    return basis


def matmul_gf2(left: List[List[int]], right: List[List[int]]) -> List[List[int]]:
    if not left or not right:
        return []
    n_inner = len(left[0])
    if n_inner != len(right):
        raise ValueError('GF(2) matrix dimensions do not align')
    out_cols = len(right[0])
    result = [[0] * out_cols for _ in range(len(left))]
    for row_idx, row in enumerate(left):
        one_positions = [idx for idx, bit in enumerate(row) if bit]
        for out_col in range(out_cols):
            accum = 0
            for pos in one_positions:
                accum ^= right[pos][out_col]
            result[row_idx][out_col] = accum
    return result


def row_reduce_masks(rows: List[int], nbits: int) -> List[int]:
    reduced = [row for row in rows if row]
    pivot_row = 0
    for col in range(nbits - 1, -1, -1):
        pivot = None
        for idx in range(pivot_row, len(reduced)):
            if (reduced[idx] >> col) & 1:
                pivot = idx
                break
        if pivot is None:
            continue
        reduced[pivot_row], reduced[pivot] = reduced[pivot], reduced[pivot_row]
        for idx in range(len(reduced)):
            if idx != pivot_row and ((reduced[idx] >> col) & 1):
                reduced[idx] ^= reduced[pivot_row]
        pivot_row += 1
    return [row for row in reduced if row]


def matrix_to_masks(matrix: List[List[int]]) -> List[int]:
    if not matrix:
        return []
    masks: List[int] = []
    for col in range(len(matrix[0])):
        mask = 0
        for row in range(len(matrix)):
            if matrix[row][col]:
                mask |= 1 << row
        masks.append(mask)
    return masks


def masks_to_matrix(masks: List[int], nbits: int) -> List[List[int]]:
    matrix = [[0] * len(masks) for _ in range(nbits)]
    for col, mask in enumerate(masks):
        for row in range(nbits):
            if (mask >> row) & 1:
                matrix[row][col] = 1
    return matrix


def compact_matrix_columns(matrix: List[List[int]], nbits: int) -> List[List[int]]:
    if not matrix:
        return matrix
    compact_masks = row_reduce_masks(matrix_to_masks(matrix), nbits)
    compact_masks.sort(key=lambda mask: (mask.bit_count(), mask))
    return masks_to_matrix(compact_masks, nbits) if compact_masks else []


def apply_matrix_row(row_bits: List[int], matrix: List[List[int]]) -> Tuple[int, ...]:
    if not matrix:
        return tuple()
    code = [0] * len(matrix[0])
    one_positions = [idx for idx, bit in enumerate(row_bits) if bit]
    for col in range(len(matrix[0])):
        accum = 0
        for pos in one_positions:
            accum ^= matrix[pos][col]
        code[col] = accum
    return tuple(code)


def compute_bank_codes(representatives: List[List[int]], matrix: List[List[int]]) -> List[Tuple[int, ...]]:
    return [apply_matrix_row(rep, matrix) for rep in representatives]


def verify_addresses(bank_sets: Sequence[BankSet], matrix: List[List[int]], bank_codes: List[Tuple[int, ...]]) -> Tuple[int, int]:
    mismatches = 0
    total = 0
    for bank_idx, bank in enumerate(bank_sets):
        expected = bank_codes[bank_idx]
        for row_bits in bank.rows:
            total += 1
            if apply_matrix_row(row_bits, matrix) != expected:
                mismatches += 1
    return mismatches, total


def remove_constant_bank_bits(matrix: List[List[int]], bank_sets: Sequence[BankSet]) -> List[List[int]]:
    if not matrix:
        return matrix
    total_addresses = sum(len(bank.addrs) for bank in bank_sets)
    keep_cols: List[int] = []
    for col in range(len(matrix[0])):
        ones = 0
        single_col = [[row[col]] for row in matrix]
        for bank in bank_sets:
            for row_bits in bank.rows:
                ones += apply_matrix_row(row_bits, single_col)[0]
        if 0 < ones < total_addresses:
            keep_cols.append(col)
    if len(keep_cols) == len(matrix[0]):
        return matrix
    return [[row[col] for col in keep_cols] for row in matrix]


def select_unique_columns(z_matrix: List[List[int]], need_bits: int) -> List[int]:
    if not z_matrix or not z_matrix[0]:
        return []
    picked: List[int] = []
    max_cols = len(z_matrix[0])
    while len(picked) < max_cols:
        best_choice = None
        best_score = None
        for col in range(max_cols):
            if col in picked:
                continue
            trial_cols = picked + [col]
            trial = [[row[idx] for idx in trial_cols] for row in z_matrix]
            rank = mat_rank_gf2(trial)
            unique_codes = len({tuple(row) for row in trial})
            score = (unique_codes, rank, -len(trial_cols), -col)
            if best_score is None or score > best_score:
                best_score = score
                best_choice = col
        if best_choice is None:
            break
        picked.append(best_choice)
        unique_codes = len({tuple(row[idx] for idx in picked) for row in z_matrix})
        if len(picked) >= need_bits and unique_codes == len(z_matrix):
            break
    return picked


def solve_region_progress(bank_sets: Sequence[BankSet], lowbit: int, highbit: int) -> RegionProgress:
    nbits = highbit - lowbit + 1
    diff_matrix: List[List[int]] = []
    representatives: List[List[int]] = []
    for bank in bank_sets:
        representatives.append(bank.rows[0][:])
        base_row = bank.rows[0]
        for row_bits in bank.rows[1:]:
            diff = [row_bits[idx] ^ base_row[idx] for idx in range(nbits)]
            if any(diff):
                diff_matrix.append(diff)

    basis_columns = nullspace_basis_gf2(diff_matrix, nbits)
    if not basis_columns:
        raise RuntimeError('No nullspace remains; widen the PA bit window.')

    nmat = [[basis_columns[col][row] for col in range(len(basis_columns))] for row in range(nbits)]
    z_matrix = matmul_gf2(representatives, nmat)
    need_bits = max(math.ceil(math.log2(max(len(bank_sets), 1))), 1)
    picked_cols = select_unique_columns(z_matrix, need_bits)
    if not picked_cols:
        picked_cols = list(range(min(len(z_matrix[0]), need_bits)))

    while len({tuple(row[idx] for idx in picked_cols) for row in z_matrix}) < len(bank_sets) and len(picked_cols) < len(z_matrix[0]):
        for col in range(len(z_matrix[0])):
            if col not in picked_cols:
                picked_cols.append(col)
                break

    matrix = [[nmat[row][col] for col in picked_cols] for row in range(nbits)]
    matrix = remove_constant_bank_bits(matrix, bank_sets)
    matrix = compact_matrix_columns(matrix, nbits)
    bank_codes = compute_bank_codes(representatives, matrix)
    unique_bank_codes = len(set(bank_codes))
    mismatches, total = verify_addresses(bank_sets, matrix, bank_codes)
    result = RegionSolveResult(
        selector=None,
        bank_sets=list(bank_sets),
        matrix=matrix,
        bank_codes=bank_codes,
        mismatches=mismatches,
        total_addresses=total,
    )
    return RegionProgress(
        result=result,
        unique_bank_codes=unique_bank_codes,
        unresolved_sets=len(bank_sets) - unique_bank_codes,
    )


def solve_region(bank_sets: Sequence[BankSet], lowbit: int, highbit: int) -> RegionSolveResult:
    progress = solve_region_progress(bank_sets, lowbit, highbit)
    if progress.unique_bank_codes != len(bank_sets):
        raise RuntimeError(
            f'Only {progress.unique_bank_codes}/{len(bank_sets)} unique bank codes were recovered; '
            'the region still contains multiple local mappings.'
        )
    return progress.result


def parse_region_file(path: str) -> List[RegionSelector]:
    selectors: List[RegionSelector] = []
    current_terms: Optional[List[SelectorTerm]] = None
    with open(path, 'r', encoding='utf-8') as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if parts[0] == 'region':
                if current_terms is not None:
                    selectors.append(RegionSelector(tuple(current_terms)))
                if len(parts) == 1:
                    current_terms = []
                elif len(parts) == 3:
                    selectors.append(RegionSelector.from_mask_value(int(parts[1], 0), int(parts[2], 0)))
                    current_terms = None
                else:
                    raise ValueError(f'{path}:{lineno}: expected "region" or "region <mask> <value>"')
                continue
            if parts[0] == 'selector':
                if current_terms is None:
                    current_terms = []
                if len(parts) < 3:
                    raise ValueError(f'{path}:{lineno}: expected "selector <0|1> <bit> [<bit> ...]"')
                expected = int(parts[1], 0)
                bits = tuple(int(part, 0) for part in parts[2:])
                current_terms.append(SelectorTerm(expected=expected, bits=bits))
                continue
            raise ValueError(f'{path}:{lineno}: unexpected line: {line}')
    if current_terms is not None:
        selectors.append(RegionSelector(tuple(current_terms)))
    if not selectors:
        raise ValueError(f'No region selectors found in {path}')
    return selectors


def assign_manual_regions(bank_sets: Sequence[BankSet], selectors: Sequence[RegionSelector]) -> List[Tuple[RegionSelector, List[BankSet]]]:
    groups: Dict[RegionSelector, List[BankSet]] = {selector: [] for selector in selectors}
    for bank in bank_sets:
        matches = [selector for selector in selectors if all(selector.matches(addr) for addr in bank.addrs)]
        if len(matches) != 1:
            raise ValueError(
                f'Bank set {bank.name} matched {len(matches)} regions; '
                'manual selectors must partition the discovered sets exactly.'
            )
        groups[matches[0]].append(bank)
    return [(selector, groups[selector]) for selector in selectors if groups[selector]]


def solve_partition(grouped_bank_sets: Sequence[Tuple[RegionSelector, List[BankSet]]], lowbit: int, highbit: int) -> List[RegionSolveResult]:
    results: List[RegionSolveResult] = []
    for selector, bank_group in grouped_bank_sets:
        solved = solve_region(bank_group, lowbit, highbit)
        solved.selector = selector
        results.append(solved)
    return results


def candidate_score(results: Sequence[RegionSolveResult]) -> Tuple[int, int, int, int]:
    total_mismatches = sum(result.mismatches for result in results)
    max_selector_complexity = max((result.selector.complexity() if result.selector is not None else 0) for result in results)
    total_selector_complexity = sum((result.selector.complexity() if result.selector is not None else 0) for result in results)
    smallest_region = min((len(result.bank_sets) for result in results), default=0)
    return (total_mismatches, max_selector_complexity, total_selector_complexity, -smallest_region)


def format_selector_bits(bits: Tuple[int, ...]) -> str:
    return ' ^ '.join(f'a{bit}' for bit in bits)


def build_branch_diagnostic(
    selector: RegionSelector,
    bank_sets: Sequence[BankSet],
    lowbit: int,
    highbit: int,
    max_term_width: int,
    used_terms: Tuple[Tuple[int, ...], ...],
    reason: str,
    diagnostics_limit: int,
) -> BranchDiagnostic:
    progress = solve_region_progress(bank_sets, lowbit, highbit)
    suggested: List[DiagnosticSplit] = []
    for bits, groups in stable_selector_terms(bank_sets, lowbit, highbit, max_term_width, used_terms):
        child_progress = [solve_region_progress(groups[idx], lowbit, highbit) for idx in (0, 1)]
        suggested.append(
            DiagnosticSplit(
                bits=bits,
                child_sizes=(len(groups[0]), len(groups[1])),
                child_uniques=(child_progress[0].unique_bank_codes, child_progress[1].unique_bank_codes),
            )
        )
        if len(suggested) >= diagnostics_limit:
            break
    return BranchDiagnostic(
        selector=selector,
        bank_set_names=[bank.name for bank in bank_sets],
        total_sets=len(bank_sets),
        unique_bank_codes=progress.unique_bank_codes,
        reason=reason,
        suggested_splits=suggested,
    )


def stable_selector_terms(
    bank_sets: Sequence[BankSet],
    lowbit: int,
    highbit: int,
    max_term_width: int,
    used_terms: Tuple[Tuple[int, ...], ...],
) -> List[Tuple[Tuple[int, ...], Dict[int, List[BankSet]]]]:
    used = set(used_terms)
    candidates: List[Tuple[Tuple[int, ...], Dict[int, List[BankSet]]]] = []
    search_bits = list(range(lowbit, highbit + 1))
    for width in range(1, max_term_width + 1):
        for bits in itertools.combinations(search_bits, width):
            if bits in used:
                continue
            groups: Dict[int, List[BankSet]] = {0: [], 1: []}
            valid = True
            for bank in bank_sets:
                values = {sum((addr >> bit) & 1 for bit in bits) & 1 for addr in bank.addrs}
                if len(values) != 1:
                    valid = False
                    break
                groups[next(iter(values))].append(bank)
            if valid and groups[0] and groups[1]:
                candidates.append((bits, groups))
    candidates.sort(key=lambda entry: (len(entry[0]), entry[0]))
    return candidates


def extend_selector(selector: RegionSelector, bits: Tuple[int, ...], expected: int) -> RegionSelector:
    return RegionSelector(selector.terms + (SelectorTerm(expected=expected, bits=bits),))


def score_split_candidate(
    groups: Dict[int, List[BankSet]],
    lowbit: int,
    highbit: int,
) -> Tuple[int, int, int, int]:
    progress = [solve_region_progress(groups[idx], lowbit, highbit) for idx in (0, 1)]
    unresolved_total = sum(item.unresolved_sets for item in progress)
    worst_child = max(item.unresolved_sets for item in progress)
    imbalance = abs(len(groups[0]) - len(groups[1]))
    solved_children = sum(1 for idx, item in enumerate(progress) if item.unique_bank_codes == len(groups[idx]))
    return (unresolved_total, worst_child, -solved_children, imbalance)


def merge_diagnostics(*diagnostic_lists: List[BranchDiagnostic]) -> List[BranchDiagnostic]:
    merged: List[BranchDiagnostic] = []
    for items in diagnostic_lists:
        merged.extend(items)
    merged.sort(key=lambda item: (item.total_sets - item.unique_bank_codes, -item.total_sets, item.selector.complexity()))
    return merged


def print_branch_diagnostics(diagnostics: Sequence[BranchDiagnostic], limit: int) -> None:
    if not diagnostics:
        return
    print('=== Unresolved Auto-Region Branches ===', file=sys.stderr)
    for diag in diagnostics[:limit]:
        print(
            f'- {diag.selector.describe()}: {diag.unique_bank_codes}/{diag.total_sets} unique bank codes; {diag.reason}',
            file=sys.stderr,
        )
        print(f'  sets: {", ".join(diag.bank_set_names)}', file=sys.stderr)
        if diag.suggested_splits:
            for split in diag.suggested_splits:
                print(
                    f'  try {format_selector_bits(split.bits)} -> '
                    f'sizes {split.child_sizes[0]}/{split.child_sizes[1]}, '
                    f'uniques {split.child_uniques[0]}/{split.child_sizes[0]} and '
                    f'{split.child_uniques[1]}/{split.child_sizes[1]}',
                    file=sys.stderr,
                )


def print_auto_trace(verbose: bool, depth: int, trace_depth_limit: int, message: str) -> None:
    if not verbose or depth > trace_depth_limit:
        return
    indent = '  ' * depth
    print(f'{indent}{message}', file=sys.stderr)


def find_auto_regions_recursive(
    bank_sets: Sequence[BankSet],
    lowbit: int,
    highbit: int,
    max_term_width: int,
    remaining_depth: int,
    selector: RegionSelector,
    max_top_splits: int,
    diagnostics_limit: int,
    verbose: bool,
    trace_depth_limit: int,
    depth: int,
) -> Optional[AutoRegionCandidate]:
    progress = solve_region_progress(bank_sets, lowbit, highbit)
    if progress.unique_bank_codes == len(bank_sets):
        solved = progress.result
        solved.selector = selector
        print_auto_trace(
            verbose,
            depth,
            trace_depth_limit,
            f'solved {selector.describe()}: {len(bank_sets)} set(s), {solved.bit_count} local bit(s)',
        )
        return AutoRegionCandidate(results=[solved], score=candidate_score([solved]), diagnostics=[])

    reason = (
        f'Only {progress.unique_bank_codes}/{len(bank_sets)} unique bank codes were recovered; '
        'the region still contains multiple local mappings.'
    )
    if remaining_depth <= 0:
        diagnostic = build_branch_diagnostic(
            selector,
            bank_sets,
            lowbit,
            highbit,
            max_term_width,
            tuple(term.bits for term in selector.terms),
            reason,
            diagnostics_limit,
        )
        print_auto_trace(verbose, depth, trace_depth_limit, f'stopped {selector.describe()}: {reason}')
        return AutoRegionCandidate(results=[], score=(sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize), diagnostics=[diagnostic])

    used_terms = tuple(term.bits for term in selector.terms)
    best: Optional[AutoRegionCandidate] = None
    split_candidates = stable_selector_terms(bank_sets, lowbit, highbit, max_term_width, used_terms)
    split_candidates.sort(key=lambda item: score_split_candidate(item[1], lowbit, highbit) + (len(item[0]), item[0]))
    split_candidates = split_candidates[:max_top_splits]
    if not split_candidates:
        diagnostic = build_branch_diagnostic(
            selector,
            bank_sets,
            lowbit,
            highbit,
            max_term_width,
            used_terms,
            'no stable selector terms remained',
            diagnostics_limit,
        )
        print_auto_trace(verbose, depth, trace_depth_limit, f'no further stable splits under {selector.describe()}')
        return AutoRegionCandidate(results=[], score=(sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize), diagnostics=[diagnostic])

    print_auto_trace(
        verbose,
        depth,
        trace_depth_limit,
        f'search {selector.describe()}: {progress.unique_bank_codes}/{len(bank_sets)} unique codes, trying {len(split_candidates)} split(s)',
    )

    for bits, groups in split_candidates:
        split_score = score_split_candidate(groups, lowbit, highbit)
        print_auto_trace(
            verbose,
            depth,
            trace_depth_limit,
            f'try {format_selector_bits(bits)} -> sizes {len(groups[0])}/{len(groups[1])}, unresolved score {split_score[0]}/{split_score[1]}',
        )
        child_results: List[RegionSolveResult] = []
        child_diags: List[BranchDiagnostic] = []
        for expected in (0, 1):
            child_selector = extend_selector(selector, bits, expected)
            child_candidate = find_auto_regions_recursive(
                groups[expected],
                lowbit,
                highbit,
                max_term_width,
                remaining_depth - 1,
                child_selector,
                max_top_splits,
                diagnostics_limit,
                verbose,
                trace_depth_limit,
                depth + 1,
            )
            if child_candidate is None:
                child_results = []
                break
            child_results.extend(child_candidate.results)
            child_diags = merge_diagnostics(child_diags, child_candidate.diagnostics)
        if not child_results:
            continue
        candidate = AutoRegionCandidate(
            results=child_results,
            score=candidate_score(child_results),
            diagnostics=child_diags,
        )
        print_auto_trace(verbose, depth, trace_depth_limit, f'candidate {format_selector_bits(bits)} score={candidate.score}')
        if best is None or candidate.score < best.score:
            best = candidate
    if best is not None:
        return best

    diagnostic = build_branch_diagnostic(
        selector,
        bank_sets,
        lowbit,
        highbit,
        max_term_width,
        used_terms,
        'all candidate splits left at least one unresolved child branch',
        diagnostics_limit,
    )
    return AutoRegionCandidate(results=[], score=(sys.maxsize, sys.maxsize, sys.maxsize, sys.maxsize), diagnostics=[diagnostic])


def find_auto_regions(
    bank_sets: Sequence[BankSet],
    lowbit: int,
    highbit: int,
    max_term_width: int,
    max_depth: int,
    max_regions: int,
    max_top_splits: int,
    diagnostics_limit: int,
    verbose: bool,
    trace_depth_limit: int,
) -> Optional[AutoRegionCandidate]:
    candidate = find_auto_regions_recursive(
        bank_sets,
        lowbit,
        highbit,
        max_term_width,
        max_depth,
        RegionSelector(tuple()),
        max_top_splits,
        diagnostics_limit,
        verbose,
        trace_depth_limit,
        0,
    )
    if candidate is None or len(candidate.results) > max_regions:
        return None
    return candidate


def format_function_lines(matrix: List[List[int]], lowbit: int) -> List[str]:
    if not matrix:
        return []
    lines: List[str] = []
    for col in range(len(matrix[0])):
        bits = [str(lowbit + row) for row in range(len(matrix)) if matrix[row][col]]
        if bits:
            lines.append(' '.join(bits))
    return lines


def write_mapping_file(path: str, results: Sequence[RegionSolveResult], lowbit: int) -> None:
    max_local_bits = max((result.bit_count for result in results), default=0)
    with open(path, 'w', encoding='utf-8') as handle:
        if len(results) == 1 and (results[0].selector is None or not results[0].selector.terms):
            for line in format_function_lines(results[0].matrix, lowbit):
                handle.write(f'{line}\n')
            return

        handle.write('# Piecewise DRAM bank mapping\n')
        handle.write(f'# global_color = (region_index << {max_local_bits}) | local_color\n')
        for idx, result in enumerate(results):
            selector = result.selector if result.selector is not None else RegionSelector(tuple())
            for line in selector.format_block():
                handle.write(f'{line}\n')
            handle.write(f'# region_index {idx}\n')
            for line in format_function_lines(result.matrix, lowbit):
                handle.write(f'{line}\n')
            handle.write('\n')


def write_regions_file(path: str, results: Sequence[RegionSolveResult]) -> None:
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('# Starter manual region file generated by gf2_bank_solver.py\n')
        handle.write('# Edit selector blocks as needed and pass this file back with --regions-file.\n\n')
        for idx, result in enumerate(results):
            selector = result.selector if result.selector is not None else RegionSelector(tuple())
            for line in selector.format_block():
                handle.write(f'{line}\n')
            if idx + 1 != len(results):
                handle.write('\n')


def print_solution_report(results: Sequence[RegionSolveResult], lowbit: int, source: str) -> None:
    print(f'=== Solve mode: {source} ===')
    max_local_bits = max((result.bit_count for result in results), default=0)
    for idx, result in enumerate(results):
        label = result.selector.describe() if result.selector is not None else 'global'
        print(f'Region {idx}: {label}')
        print(f'  bank sets: {len(result.bank_sets)}')
        print(f'  local bank bits: {result.bit_count}')
        print(f'  verification: {result.total_addresses - result.mismatches}/{result.total_addresses} correct')
        if len(results) > 1:
            print(f'  global color encoding: (region_index << {max_local_bits}) | local_color')
        function_lines = format_function_lines(result.matrix, lowbit)
        for col, line in enumerate(function_lines):
            pretty = ' ^ '.join(f'a{bit}' for bit in map(int, line.split()))
            print(f'  bank_bit[{col}] = {pretty}')
        if not function_lines:
            print('  bank_bit[*] = <none>')
        print('')


def main() -> None:
    parser = argparse.ArgumentParser(description='GF(2) DRAM bank mapping solver from per-bank address files.')
    parser.add_argument('--files', nargs='+', required=True, help='One file per bank (each contains PAs for that bank).')
    parser.add_argument('--lowbit', type=int, default=5, help='Lowest PA bit to include (default: 5).')
    parser.add_argument('--highbit', type=int, default=40, help='Highest PA bit to include (inclusive).')
    parser.add_argument('--output', default='recovered_bank_mapping.txt', help='Output mapping file path.')
    parser.add_argument('--regions-file', help='Manual region definitions using region/selector blocks or region <mask> <value>.')
    parser.add_argument('--emit-regions-file', help='Write the recovered selectors as a starter manual regions file.')
    parser.add_argument('--auto-regions', action='store_true', help='Infer a piecewise mapping by recursively searching for selector terms.')
    parser.add_argument('--auto-max-selector-bits', type=int, default=2, help='Maximum number of PA bits in one auto-discovered selector term.')
    parser.add_argument('--auto-max-depth', type=int, default=4, help='Maximum recursive auto-region split depth.')
    parser.add_argument('--auto-max-regions', type=int, default=8, help='Maximum number of regions considered during auto-discovery.')
    parser.add_argument('--auto-top-splits', type=int, default=16, help='Maximum number of promising selector terms explored per unresolved branch.')
    parser.add_argument('--auto-diagnostics-limit', type=int, default=8, help='Maximum number of unresolved branches or suggested splits reported on failure.')
    parser.add_argument('--auto-trace-depth', type=int, default=0, help='Maximum recursion depth printed by --verbose auto-region tracing.')
    parser.add_argument('--verbose', action='store_true', help='Print extra diagnostics.')
    args = parser.parse_args()

    if args.highbit < args.lowbit:
        print('highbit must be >= lowbit', file=sys.stderr)
        sys.exit(1)

    bank_sets = load_bank_sets(args.files, args.lowbit, args.highbit)
    if args.verbose:
        print(
            f'Loaded {len(bank_sets)} bank sets using PA bits a{args.lowbit}..a{args.highbit}',
            file=sys.stderr,
        )

    source = 'global'
    if args.regions_file:
        selectors = parse_region_file(args.regions_file)
        results = solve_partition(assign_manual_regions(bank_sets, selectors), args.lowbit, args.highbit)
        source = 'manual-regions'
    elif args.auto_regions:
        candidate = find_auto_regions(
            bank_sets,
            args.lowbit,
            args.highbit,
            args.auto_max_selector_bits,
            args.auto_max_depth,
            args.auto_max_regions,
            args.auto_top_splits,
            args.auto_diagnostics_limit,
            args.verbose,
            args.auto_trace_depth,
        )
        if candidate is not None and candidate.results and candidate.score[0] == 0:
            results = candidate.results
            source = 'auto-regions'
        else:
            if candidate is not None:
                print_branch_diagnostics(candidate.diagnostics, args.auto_diagnostics_limit)
            try:
                results = [solve_region(bank_sets, args.lowbit, args.highbit)]
            except RuntimeError as exc:
                print(f'Auto-region search failed and no unique global mapping exists: {exc}', file=sys.stderr)
                sys.exit(2)
            if candidate is not None and args.verbose:
                print(f'Falling back to global solve; best auto score was {candidate.score}', file=sys.stderr)
    else:
        try:
            results = [solve_region(bank_sets, args.lowbit, args.highbit)]
        except RuntimeError as exc:
            print(f'Global solve failed: {exc}', file=sys.stderr)
            sys.exit(2)

    if len(results) == 1 and source == 'global':
        results[0].selector = None

    print_solution_report(results, args.lowbit, source)
    write_mapping_file(args.output, results, args.lowbit)
    print(f'Stored recovered mapping to {args.output}')
    if args.emit_regions_file:
        write_regions_file(args.emit_regions_file, results)
        print(f'Stored starter regions file to {args.emit_regions_file}')


if __name__ == '__main__':
    main()
