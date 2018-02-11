#include "catch.hpp"
#include <squanch/Model.h>
#include <squanch/Line.h>

TEST_CASE("Model", "[model]")
{
  SECTION("clear") 
  {
    using namespace squanch;
    using namespace Eigen;

    squanch::Model<double> model;

    Line<double, Vector3d> line;
    line.p1 = Vector3d(0, 0, 0);
    line.p2 = Vector3d(10, 0, 0);

    auto&& nline1 = line.to_nurbs(model);
    auto&& nline2 = line.to_nurbs(model);
    auto&& nline3 = line.to_nurbs(model);

    REQUIRE(model.patches().size() == 3);
    REQUIRE(model.num_nodes() == 6);

    model.clear();

    REQUIRE(model.patches().size() == 0);
    REQUIRE(model.num_nodes() == 0);
  }
}
