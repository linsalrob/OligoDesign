"""Core DNA sequence object for the OligoDesign library."""

from __future__ import annotations

import math
from typing import Iterator

# Standard unambiguous DNA bases
VALID_BASES: frozenset[str] = frozenset("ACGT")

# IUPAC ambiguity codes (includes the four standard bases)
IUPAC_BASES: frozenset[str] = frozenset("ACGTRYMKSWHBVDN")

# Translation table for complementing uppercase DNA (ACGT only)
_COMPLEMENT_TABLE: dict[int, int] = str.maketrans("ACGT", "TGCA")

# IUPAC complement translation (covers standard bases and all ambiguity codes)
_IUPAC_COMPLEMENT_TABLE: dict[int, int] = str.maketrans(
    "ACGTRYMKSWHBVDN",
    "TGCAYRKMSWDVBHN",
)


class DNA:
    """A DNA sequence object.

    Stores a DNA sequence in uppercase and provides common sequence
    operations including reverse complementing, alphabet validation,
    GC content, motif search, and structural heuristics.

    Indexing follows standard Python 0-based conventions. For 1-based
    indexing use the explicit methods :meth:`get_base_1` and
    :meth:`slice_1`.

    Parameters
    ----------
    sequence:
        The DNA sequence string. Case-insensitive; stored as uppercase.
    allow_ambiguous:
        If ``True``, accept IUPAC ambiguity codes beyond ``ACGT``.
        Default is ``False``.

    Raises
    ------
    TypeError
        If *sequence* is not a string.
    ValueError
        If *sequence* contains characters outside the accepted alphabet.

    Examples
    --------
    >>> dna = DNA("ACGT")
    >>> str(dna)
    'ACGT'
    >>> str(dna.reverse_complement())
    'ACGT'
    >>> dna.gc_content()
    0.5
    """

    def __init__(self, sequence: str, allow_ambiguous: bool = False) -> None:
        if not isinstance(sequence, str):
            raise TypeError(
                f"DNA sequence must be a string, not {type(sequence).__name__!r}"
            )

        upper = sequence.upper()
        valid = IUPAC_BASES if allow_ambiguous else VALID_BASES
        invalid = frozenset(upper) - valid
        if invalid:
            sorted_invalid = sorted(invalid)
            alphabet_name = "IUPAC bases" if allow_ambiguous else "ACGT"
            raise ValueError(
                f"Invalid DNA characters: {sorted_invalid!r}. "
                f"Expected {alphabet_name}."
            )

        self._sequence: str = upper
        self._allow_ambiguous: bool = allow_ambiguous

    # ------------------------------------------------------------------
    # Python data model
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._sequence)

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        return f"DNA({self._sequence!r})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, DNA):
            return self._sequence == other._sequence
        if isinstance(other, str):
            return self._sequence == other.upper()
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._sequence)

    def __iter__(self) -> Iterator[str]:
        return iter(self._sequence)

    def __getitem__(self, key: int | slice) -> str:
        """Return base(s) at the given 0-based index or slice."""
        return self._sequence[key]

    def __contains__(self, item: object) -> bool:
        if isinstance(item, DNA):
            return str(item) in self._sequence
        if isinstance(item, str):
            return item.upper() in self._sequence
        return False

    # ------------------------------------------------------------------
    # 1-based indexing
    # ------------------------------------------------------------------

    def get_base_1(self, position: int) -> str:
        """Return the base at the given 1-based position.

        Parameters
        ----------
        position:
            1-based position (1 = first base).

        Returns
        -------
        str
            Single character base at *position*.

        Raises
        ------
        IndexError
            If *position* is outside ``[1, len(self)]``.

        Examples
        --------
        >>> DNA("ACGT").get_base_1(1)
        'A'
        >>> DNA("ACGT").get_base_1(4)
        'T'
        """
        n = len(self._sequence)
        if position < 1 or position > n:
            raise IndexError(
                f"Position {position} is out of range for a sequence of length {n}. "
                f"Valid 1-based positions are 1 to {n}."
            )
        return self._sequence[position - 1]

    def slice_1(self, start: int, end: int) -> DNA:
        """Return a subsequence using 1-based, inclusive coordinates.

        Parameters
        ----------
        start:
            Start position (1-based, inclusive).
        end:
            End position (1-based, inclusive).

        Returns
        -------
        DNA
            New :class:`DNA` object containing the subsequence.

        Raises
        ------
        IndexError
            If *start* or *end* are outside the valid range.
        ValueError
            If *start* > *end*.

        Examples
        --------
        >>> str(DNA("ACGTACGT").slice_1(1, 4))
        'ACGT'
        >>> str(DNA("ACGTACGT").slice_1(5, 8))
        'ACGT'
        """
        n = len(self._sequence)
        if start < 1 or end > n:
            raise IndexError(
                f"Slice [{start}:{end}] is out of range for a sequence of length {n}. "
                f"Valid 1-based positions are 1 to {n}."
            )
        if start > end:
            raise ValueError(
                f"Start position {start} must not be greater than end position {end}."
            )
        return DNA(self._sequence[start - 1 : end], allow_ambiguous=self._allow_ambiguous)

    # ------------------------------------------------------------------
    # Sequence operations
    # ------------------------------------------------------------------

    def reverse(self) -> DNA:
        """Return the reverse of the sequence (5'->3' order reversed).

        Returns
        -------
        DNA
            New :class:`DNA` object with the sequence reversed.

        Examples
        --------
        >>> str(DNA("ACGT").reverse())
        'TGCA'
        """
        return DNA(self._sequence[::-1], allow_ambiguous=self._allow_ambiguous)

    def complement(self) -> DNA:
        """Return the complement of the sequence (5'->3' direction unchanged).

        Each base is replaced by its Watson-Crick complement. IUPAC
        ambiguity codes are complemented correctly when
        ``allow_ambiguous`` was set.

        Returns
        -------
        DNA
            New :class:`DNA` object with complemented bases.

        Examples
        --------
        >>> str(DNA("ACGT").complement())
        'TGCA'
        """
        table = _IUPAC_COMPLEMENT_TABLE if self._allow_ambiguous else _COMPLEMENT_TABLE
        return DNA(self._sequence.translate(table), allow_ambiguous=self._allow_ambiguous)

    def reverse_complement(self) -> DNA:
        """Return the reverse complement of the sequence.

        Equivalent to reading the opposite strand 5'->3'.

        Returns
        -------
        DNA
            New :class:`DNA` object containing the reverse complement.

        Examples
        --------
        >>> str(DNA("GAATTC").reverse_complement())
        'GAATTC'
        >>> str(DNA("AAAA").reverse_complement())
        'TTTT'
        """
        table = _IUPAC_COMPLEMENT_TABLE if self._allow_ambiguous else _COMPLEMENT_TABLE
        return DNA(
            self._sequence[::-1].translate(table),
            allow_ambiguous=self._allow_ambiguous,
        )

    # ------------------------------------------------------------------
    # Composition and content
    # ------------------------------------------------------------------

    def gc_content(self) -> float:
        """Return GC content as a fraction between 0.0 and 1.0.

        Returns 0.0 for an empty sequence.

        Returns
        -------
        float
            Fraction of G and C bases in the sequence.

        Examples
        --------
        >>> DNA("ACGT").gc_content()
        0.5
        >>> DNA("GCGC").gc_content()
        1.0
        """
        if not self._sequence:
            return 0.0
        gc = self._sequence.count("G") + self._sequence.count("C")
        return gc / len(self._sequence)

    def base_composition(self) -> dict[str, int]:
        """Return a dictionary of base counts for ACGT.

        Returns
        -------
        dict[str, int]
            Mapping of ``{base: count}`` for each of A, C, G, T.

        Examples
        --------
        >>> DNA("AACGTT").base_composition()
        {'A': 2, 'C': 1, 'G': 1, 'T': 2}
        """
        return {base: self._sequence.count(base) for base in "ACGT"}

    def entropy(self) -> float:
        """Return the Shannon entropy of the sequence in bits.

        Computed over the four canonical bases (A, C, G, T).  Returns 0.0
        for an empty sequence or a single-base sequence.

        The maximum value is ``log2(4) = 2.0`` bits, achieved when all four
        bases are equally represented.

        Returns
        -------
        float
            Shannon entropy in bits.

        Examples
        --------
        >>> DNA("ACGT").entropy()
        2.0
        >>> DNA("AAAA").entropy()
        0.0
        """
        n = len(self._sequence)
        if n == 0:
            return 0.0
        result = 0.0
        for base in "ACGT":
            count = self._sequence.count(base)
            if count > 0:
                p = count / n
                result -= p * math.log2(p)
        return result

    # ------------------------------------------------------------------
    # Motif search
    # ------------------------------------------------------------------

    def find_motif(self, motif: str) -> list[int]:
        """Find all 0-based positions where *motif* occurs.

        Overlapping matches are reported. The search is case-insensitive.

        Parameters
        ----------
        motif:
            DNA motif to search for.

        Returns
        -------
        list[int]
            Sorted list of 0-based start positions where *motif* begins.

        Examples
        --------
        >>> DNA("ACGTACGT").find_motif("ACG")
        [0, 4]
        >>> DNA("AAAAAA").find_motif("AA")
        [0, 1, 2, 3, 4]
        """
        motif_upper = motif.upper()
        positions: list[int] = []
        start = 0
        while True:
            pos = self._sequence.find(motif_upper, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        return positions

    # ------------------------------------------------------------------
    # Sequence quality checks
    # ------------------------------------------------------------------

    def longest_homopolymer(self) -> int:
        """Return the length of the longest homopolymer run.

        A homopolymer is a consecutive run of identical nucleotides
        (e.g. ``AAAA``).

        Returns
        -------
        int
            Length of the longest run, or 0 for an empty sequence.

        Examples
        --------
        >>> DNA("ACAAAAT").longest_homopolymer()
        4
        >>> DNA("ACGT").longest_homopolymer()
        1
        """
        if not self._sequence:
            return 0
        max_run = 1
        current_run = 1
        for i in range(1, len(self._sequence)):
            if self._sequence[i] == self._sequence[i - 1]:
                current_run += 1
                if current_run > max_run:
                    max_run = current_run
            else:
                current_run = 1
        return max_run

    def has_homopolymer(self, min_length: int = 4) -> bool:
        """Return ``True`` if the sequence contains a homopolymer run of at
        least *min_length* bases.

        Parameters
        ----------
        min_length:
            Minimum consecutive run length. Default is 4.

        Examples
        --------
        >>> DNA("ACAAAAT").has_homopolymer(min_length=4)
        True
        >>> DNA("ACGT").has_homopolymer(min_length=4)
        False
        """
        return self.longest_homopolymer() >= min_length

    def is_low_complexity(self, window: int = 10, threshold: float = 0.7) -> bool:
        """Return ``True`` if any window in the sequence is low complexity.

        Low complexity is defined as a single base exceeding *threshold*
        fraction of the bases within a sliding *window*.

        Parameters
        ----------
        window:
            Window size in bases. Default is 10.
        threshold:
            Dominant-base fraction at or above which a window is
            considered low complexity. Default is 0.7.

        Examples
        --------
        >>> DNA("AAAAAAAAAG").is_low_complexity(window=10, threshold=0.7)
        True
        >>> DNA("ACGTACGTAC").is_low_complexity(window=10, threshold=0.7)
        False
        """
        seq = self._sequence
        n = len(seq)
        if n == 0:
            return False
        actual_window = min(window, n)
        for i in range(n - actual_window + 1):
            sub = seq[i : i + actual_window]
            dominant = max(sub.count(b) for b in "ACGT")
            if dominant / actual_window >= threshold:
                return True
        return False

    # ------------------------------------------------------------------
    # Self-complementarity and hairpin heuristics
    # ------------------------------------------------------------------

    def is_self_complementary(self) -> bool:
        """Return ``True`` if the sequence is equal to its own reverse complement.

        A DNA sequence that equals its reverse complement is a perfect
        palindrome (e.g. the EcoRI site ``GAATTC``).

        Examples
        --------
        >>> DNA("GAATTC").is_self_complementary()
        True
        >>> DNA("ACGT").is_self_complementary()
        False
        """
        return self._sequence == str(self.reverse_complement())

    def has_hairpin(
        self,
        min_stem: int = 4,
        min_loop: int = 3,
        max_loop: int = 8,
    ) -> bool:
        """Return ``True`` if the sequence likely contains a hairpin structure.

        Uses a simple heuristic: searches for an inverted repeat separated
        by a loop of configurable length. This is **not** a thermodynamic
        calculation; it is a fast structural pre-filter.

        Parameters
        ----------
        min_stem:
            Minimum number of base pairs forming the stem. Default is 4.
        min_loop:
            Minimum loop length in bases. Default is 3.
        max_loop:
            Maximum loop length in bases. Default is 8.

        Returns
        -------
        bool
            ``True`` if an inverted repeat consistent with a hairpin stem
            of at least *min_stem* base pairs is found.

        Examples
        --------
        >>> DNA("AAAACCCTTTT").has_hairpin(min_stem=4, min_loop=3, max_loop=8)
        True
        """
        seq = self._sequence
        n = len(seq)
        for stem_len in range(min_stem, n // 2 + 1):
            for loop_len in range(min_loop, max_loop + 1):
                required = 2 * stem_len + loop_len
                if required > n:
                    break
                for i in range(n - required + 1):
                    stem5 = seq[i : i + stem_len]
                    stem3 = seq[i + stem_len + loop_len : i + required]
                    rc_stem5 = stem5[::-1].translate(_COMPLEMENT_TABLE)
                    if stem3 == rc_stem5:
                        return True
        return False
