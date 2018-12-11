#include <fstream>
#include <iostream>
#include <vector>

#include <opennurbs_public.h>

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " input.bzr output.3dm" << std::endl;
    return 1;
  }

  ON::Begin();
  ON_TextLog error_log;

  // Read the data

  int n, d, t;
  ON_BezierSurface surf;
  std::vector<ON_2dPoint> domain;
  std::vector<ON_BezierCurve> trim_curves3d;
  {
    double x, y, z, w;
    std::ifstream f(argv[1]);
    f >> n >> d >> t;
    surf.Create(3, true, t + 1, t + 1);
    for (int i = 0; i <= t; ++i)
      for (int j = 0; j <= t; ++j) {
        f >> x >> y >> z >> w;
        surf.SetCV(i, j, ON_4dPoint(x, y, z, w));
      }
    for (int i = 0; i < n; ++i) {
      f >> x >> y;
      domain.emplace_back(x, y);
    }
	for (int i = 0; i < n; ++i) {
		ON_3dPointArray cpts;
		for (int j = 0; j <= d; ++j) {
			f >> x >> y >> z;
			cpts.Append(ON_3dPoint(x, y, z));
		}
		trim_curves3d.emplace_back(cpts);
	}
    if (!f) {
      std::cerr << "Error while parsing file." << std::endl;
      return 2;
    }
  }

  // Generate 2D trim curves
  
  std::vector<ON_LineCurve *> trim_curves;
  for (int i = 0; i < n; ++i) {
    const auto &p = domain[(i+n-1)%n];
	const auto &q = domain[i];
    trim_curves.push_back(new ON_LineCurve(p, q));
  }

  // Create the B-rep topology

  ONX_Model model;
  model.m_sStartSectionComments = "This file was generated by bzr-to-3dm.";
  model.AddDefaultLayer(nullptr, ON_Color(255, 255, 255));

  ON_Brep brep;
  auto *nurbs = new ON_NurbsSurface(surf);
  int nurbs_id = brep.m_S.Count();
  brep.m_S.Append(nurbs);

  for (int i = 0; i < n; ++i) {
	  brep.m_C2.Append(trim_curves[i]);
	  brep.m_C3.Append(new ON_NurbsCurve(trim_curves3d[i]));
  }

  for (const auto &p : domain) {
    brep.NewVertex(surf.PointAt(p.x, p.y));
    brep.m_V.Last()->m_tolerance = 0.0;
  }

  for (int i = 0; i < n; ++i) {
    brep.NewEdge(brep.m_V[(i+n-1)%n], brep.m_V[i], i);
    brep.m_E.Last()->m_tolerance = 0.0;
  }

  ON_BrepFace &face = brep.NewFace(nurbs_id);
  ON_BrepLoop &loop = brep.NewLoop(ON_BrepLoop::outer, face);

  for (int i = 0; i < n; ++i) {
    brep.NewTrim(brep.m_E[i], false, loop, i);
    brep.m_T[i].m_type = ON_BrepTrim::boundary;
    brep.m_T[i].m_tolerance[0] = 0.0;
    brep.m_T[i].m_tolerance[1] = 0.0;
  }

  model.AddModelGeometryComponent(&brep, nullptr);  
  model.Write(argv[2], 0, &error_log);
}
