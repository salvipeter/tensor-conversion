#include <opennurbs_bezier.h>
#include <opennurbs_brep.h>

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " input.bzr output.3dm" << std::endl;
    return 1;
  }

  size_t n, d;
  ON_BezierSurface surf;
  ON_Curve trim;
  std::vector<ON_LineCurve> trim_curves;
  {
    double x, y, z, w;
    std::ifstream f(argv[1]);
    f >> n >> d;
    surf.Create(3, true, d + 1, d + 1);
    for (size_t i = 0; i <= d; ++i)
      for (size_t j = 0; j <= d; ++j) {
        f >> x >> y >> z >> w;
        surf.SetCV(i, j, {x, y, z, w});
      }
    std::vector<ON_2dPoint> domain;
    for (size_t i = 0; i < n; ++i) {
      f >> x >> y;
      domain.emplace_back(x, y);
    }
    for (size_t i = 0; i < n; ++i) {
      const auto &p = domain[i];
      const auto &q = domain[(i+1)%n];
      trim_curves.emplace_back(p, q);
    }
    if (!f) {
      std::cerr << "Error while parsing file." << std::endl;
      return 2;
    }
  }

  ON_Brep brep;
  ON_BrepLoop loop;
  int surf_id = brep.AddSurface(&surf);
  int trim_id = brep.AddTrimCurve(&trim);

  ONX_Model model;
  ON_3dmObjectAttributes attr;
  model.AddModelGeometryComponent(&brep, &attr);
  model.Write(argv[2]);
}
