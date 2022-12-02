/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#include "NekUtils.h"

namespace nekrs
{

std::vector<std::vector<int>>
cornerGLLIndices(const int & n, const bool exact)
{
  std::vector<std::vector<int>> corner_indices;

  if (!exact)
  {
    corner_indices.resize(1);
    int back_corner = (n + 1) * (n + 1) * n;
    corner_indices[0] = {0,
                         n,
                         (n + 1) * n,
                         (n + 1) * (n + 1) - 1,
                         back_corner,
                         back_corner + n,
                         back_corner + (n + 1) * n,
                         back_corner + (n + 1) * (n + 1) - 1};
  }
  else
  {
    for (int depth = 0; depth < n; ++depth)
    {
      for (int row = 0; row < n; ++row)
      {
        for (int col = 0; col < n; ++col)
        {
          int start = row * (n + 1) + col + depth * (n + 1) * (n + 1);
          int back_corner = (n + 1) * (n + 1) + start;
          corner_indices.push_back({start,
                                    start + 1,
                                    start + n + 1,
                                    start + n + 2,
                                    back_corner,
                                    back_corner + 1,
                                    back_corner + n + 1,
                                    back_corner + n + 2});
        }
      }
    }
  }

  return corner_indices;
}

} // end namespace nekrs
