r"""
G-Sets
"""
#*****************************************************************************
#  Copyright (C) 2008 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from .sets_cat import Sets

#############################################################
# GSets
#     `G`-Sets play an important role in permutation groups.
#############################################################
class GSets(Category):
    """
    The category of `G`-sets, for a group `G`.

    EXAMPLES::

        sage: S = SymmetricGroup(3)                                                     # optional - sage.groups
        sage: GSets(S)                                                                  # optional - sage.groups
        Category of G-sets for Symmetric group of order 3! as a permutation group

    TODO: should this derive from Category_over_base?
    """
    def __init__(self, G):
        """
        TESTS::

            sage: S8 = SymmetricGroup(8)                                                # optional - sage.groups
            sage: TestSuite(GSets(S8)).run()                                            # optional - sage.groups
        """
        Category.__init__(self)
        self.__G = G

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: GSets(SymmetricGroup(8))  # indirect doctests                         # optional - sage.groups
            Category of G-sets for Symmetric group of order 8! as a permutation group
        """
        return "G-sets for %s"%self.__G

    #def construction(self):
    #    return (self.__class__, self.__G)

    def super_categories(self):
        """
        EXAMPLES::

            sage: GSets(SymmetricGroup(8)).super_categories()                           # optional - sage.groups
            [Category of sets]
        """
        return [Sets()]

    @classmethod
    def an_instance(cls):
        """
        Returns an instance of this class.

        EXAMPLES::

            sage: GSets.an_instance()  # indirect doctest                               # optional - sage.groups
            Category of G-sets for Symmetric group of order 8! as a permutation group
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        G = SymmetricGroup(8)
        return cls(G)
