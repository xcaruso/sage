r"""
Features for testing the presence of the SageMath interfaces to ``gap`` and of GAP packages
"""
# *****************************************************************************
#       Copyright (C) 2016 Julian Rüth
#                     2018 Jeroen Demeyer
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Feature, FeatureTestResult, PythonModule
from .join_feature import JoinFeature
from .sagemath import sage__libs__gap

class GapPackage(Feature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of a GAP package.

    .. SEEALSO::

        :class:`Feature sage.libs.gap <~sage.features.sagemath.sage__libs__gap>`

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")
        Feature('gap_package_grape')
    """
    def __init__(self, package, **kwds):
        r"""
        TESTS::

            sage: from sage.features.gap import GapPackage
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), GapPackage)
            True
        """
        Feature.__init__(self, f"gap_package_{package}", **kwds)
        self.package = package

    def _is_present(self):
        r"""
        Return whether the package is available in GAP.

        This does not check whether this package is functional.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape", spkg="gap_packages")._is_present()  # optional - gap_packages
            FeatureTestResult('gap_package_grape', True)
        """
        try:
            from sage.libs.gap.libgap import libgap
        except ImportError:
            return FeatureTestResult(self, False,
                                     reason="sage.libs.gap is not available")
        command = 'TestPackageAvailability("{package}")'.format(package=self.package)
        presence = libgap.eval(command)
        if presence:
            return FeatureTestResult(self, True,
                    reason="`{command}` evaluated to `{presence}` in GAP.".format(command=command, presence=presence))
        else:
            return FeatureTestResult(self, False,
                    reason="`{command}` evaluated to `{presence}` in GAP.".format(command=command, presence=presence))


def all_features():
    return []
