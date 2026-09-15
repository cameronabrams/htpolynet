"""Bond-context descriptor classes.

These are pure data objects that describe a bond by the chemical
context (atom names, residue names, order, bystanders, one-aways)
of its two participating atoms.  They carry no cure-specific logic
and have no htpolynet dependencies.

Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
import logging

from copy import deepcopy

logger = logging.getLogger(__name__)


class BondTemplate:
    def __init__(self, names, resnames, intraresidue, order, bystander_resnames, bystander_atomnames, oneaway_resnames, oneaway_atomnames, siblings=None):
        """Creates a BondTemplate object.

        Args:
            names (list-like container of two ints): names of atoms in the bond
            resnames (list-like container of two strs): names of the two residues to which the atoms belong
            intraresidue (bool): True if this is an intraresidue bond
            order (int): bond order (1=single, 2=double, ...)
            bystander_resnames (list of two list-like containers of strs): lists of names of bystander residues (residues also bound to one of the atoms), one for each atom in names
            bystander_atomnames (list of two list-like containers of strs): lists of names of atoms in bystander residues (parallel to bystander_resnames)
            oneaway_resnames (list-like container of strs): names of one-away residues (residues bound one bond away from the new interresidue bond; only relevant for C=C free-radical polymerization)
            oneaway_atomnames (list-like container of strs): names of atoms in one-away residues
            siblings (list of two lists of (str, str) tuples): for each atom in names, the
                (atomName, partner resName) of every other atom in its residue, within two
                intraresidue bonds, that is already bonded to another residue; defaults to
                none on either side.  See TopoCoord.get_siblings.
        """
        self.names = names
        self.resnames = resnames
        self.intraresidue = intraresidue
        self.bystander_resnames = bystander_resnames
        self.bystander_atomnames = bystander_atomnames
        self.oneaway_resnames = oneaway_resnames
        self.oneaway_atomnames = oneaway_atomnames
        self.order = order
        self.siblings = [sorted(tuple(t) for t in x) for x in siblings] if siblings else [[], []]

    def reverse(self):
        """Reverses the order of all parallel lists in a BondTemplate object."""
        self.names = self.names[::-1]
        self.resnames = self.resnames[::-1]
        self.bystander_resnames = self.bystander_resnames[::-1]
        self.bystander_atomnames = self.bystander_atomnames[::-1]
        self.oneaway_resnames = self.oneaway_resnames[::-1]
        self.oneaway_atomnames = self.oneaway_atomnames[::-1]
        self.siblings = self.siblings[::-1]

    def __str__(self):
        return f'BondTemplate {self.names} resnames {self.resnames} intraresidue? {self.intraresidue} order {self.order} bystander-resnames {self.bystander_resnames} bystander-atomnames {self.bystander_atomnames} oneaway-resnames {self.oneaway_resnames} oneaway-atomnames {self.oneaway_atomnames} siblings {self.siblings}'

    def __eq__(self,other):
        check = self.names == other.names
        check = check and self.intraresidue == other.intraresidue
        check = check and self.resnames == other.resnames
        check = check and self.bystander_resnames == other.bystander_resnames
        check = check and self.bystander_atomnames == other.bystander_atomnames
        check = check and self.oneaway_resnames == other.oneaway_resnames
        check = check and self.oneaway_atomnames == other.oneaway_atomnames
        check = check and self.siblings == other.siblings
        return check

    def is_reverse_of(self,other):
        """Returns True if self and other are reverse of each other.

        Args:
            other (BondTemplate): another BondTemplate object

        Returns:
            bool: True if self and other are reverse copies of each other
        """
        rb = deepcopy(other)
        rb.reverse()
        return self == rb

    def bystander_count(self):
        """Total number of bystander (resname, atomname) pairs and siblings across both sides.

        Used by find_template's best-match selection: more context = more
        specific template, which should win when multiple templates can match
        an instance under subset semantics.  Siblings count here so that a
        triazine template with two ring carbons already bonded beats the one
        with one, which beats the bare one.
        """
        return sum(len(s) for s in self.bystander_resnames) + sum(len(s) for s in self.siblings)

    def matches(self, other):
        """Return True if `self` (a parameterization-stage template) is a
        valid description of `other` (an as-found system bond instance).

        Strict equality on:

        - atom names (``names``)
        - residue names (``resnames``) and ``intraresidue`` flag
        - one-away context (``oneaway_resnames`` / ``oneaway_atomnames``) —
          chain-extension templates use these to discriminate dimer- vs
          trimer- vs tetramer-context bonds, so the match has to be exact.

        Subset semantics on bystanders:

        - paired ``(bystander_resname, bystander_atomname)`` items are
          treated as a multiset on each side; every item the template
          declares must appear in the instance, but the instance may
          carry additional bystanders the template does not mention.

        Rationale: a small-fragment parameterization template (e.g. the
        CY-CY cyanate dimer) carries no inter-residue context at the
        cure-reactive atom because the fragment is a single residue.
        The system-instance bond, formed inside an assembled monomer
        (e.g. BCY = BPA + 2 CY), naturally picks up the BPA-O ester
        bystander at the same atom.  That extra structural context does
        not change the local chemistry of the bond template, so the
        template should still match.  Exact bystander equality would
        force the cure reactant to mirror every inter-residue partner
        the instance carries, which forecloses the small-fragment idiom.

        Subset semantics on siblings, for the same reason.  A sibling is
        another atom of a bonding residue, close enough to the bonding atom
        that the template overwrites its charge, that is already bonded to
        another residue.  A template declaring a sibling describes that atom
        as reacted, so it must not be used where the atom has not; a template
        declaring none has always been used regardless, and still is.

        Args:
            other (BondTemplate): the instance bond to test for compatibility

        Returns:
            bool: True if self can serve as a template for other
        """
        if self.names != other.names: return False
        if self.intraresidue != other.intraresidue: return False
        if self.resnames != other.resnames: return False
        if self.oneaway_resnames != other.oneaway_resnames: return False
        if self.oneaway_atomnames != other.oneaway_atomnames: return False
        for side in (0, 1):
            t_pairs = list(zip(self.bystander_resnames[side],
                               self.bystander_atomnames[side]))
            i_pairs = list(zip(other.bystander_resnames[side],
                               other.bystander_atomnames[side]))
            remaining = list(i_pairs)
            for tp in t_pairs:
                if tp in remaining:
                    remaining.remove(tp)
                else:
                    return False
            remaining = list(other.siblings[side])
            for sib in self.siblings[side]:
                if sib in remaining:
                    remaining.remove(sib)
                else:
                    return False
        return True

    def matches_reverse_of(self, other):
        """Like ``matches`` but tries the reversed orientation of ``other``."""
        rb = deepcopy(other)
        rb.reverse()
        return self.matches(rb)


BondTemplateList = list[BondTemplate]


class ReactionBond:
    def __init__(self, idx, resids, order, bystanders, bystanders_atomidx, oneaways, oneaways_atomidx):
        """Generates a new ReactionBond object.

        Args:
            idx (list-like container of two ints): indices of two atoms that form the bond
            resids (list-like container of two ints): resids of the two atoms that form the bond
            order (int): order of the bond (1=single, 2=double, ...)
            bystanders (list of two list-like containers of ints): two lists of indices of interresidue resids already bound to each atom in idx
            bystanders_atomidx (list of two list-like containers of ints): two lists of indices of interresidue atoms already bound to each atom in idx
            oneaways (list of ints): list of one-away resids
            oneaways_atomidx (list of ints): list of one-away atom indices
        """
        self.idx = idx
        self.resids = resids
        self.bystander_resids = bystanders
        self.bystander_atomidx = bystanders_atomidx
        self.oneaway_resids = oneaways
        self.oneaway_atomidx = oneaways_atomidx
        self.order = order

    def reverse(self):
        self.idx = self.idx[::-1]
        self.resids = self.resids[::-1]
        self.bystander_resids = self.bystander_resids[::-1]
        self.bystander_atomidx = self.bystander_atomidx[::-1]
        self.oneaway_resids = self.oneaway_resids[::-1]
        self.oneaway_atomidx = self.oneaway_atomidx[::-1]

    def __str__(self):
        return f'ReactionBond {self.idx} resids {self.resids} order {self.order} bystander-resids {self.bystander_resids} oneaway-resids {self.oneaway_resids} oneaway-atomidx {self.oneaway_atomidx}'


ReactionBondList = list[ReactionBond]
