#include <boost/test/unit_test.hpp>
#include <squanch/Model.h>
#include <squanch/Line.h>
#include <squanch/Extrude.h>
#include <squanch/Refinement.h>
#include <squanch/MortarConnector.h>
#include "TestUtils.h"

template <typename MatType> void checkRpMortar(const MatType& mortar)
{
	BOOST_CHECK_CLOSE(mortar(0, 0), 1, 1);
	BOOST_CHECK_SMALL(mortar(0, 1), 1e-7);
	BOOST_CHECK_SMALL(mortar(0, 2), 1e-7);
	BOOST_CHECK_SMALL(mortar(0, 3), 1e-7);
	BOOST_CHECK_SMALL(mortar(0, 4), 1e-7);
	BOOST_CHECK_SMALL(mortar(0, 5), 1e-7);

	BOOST_CHECK_CLOSE(mortar(1, 0), 0.25, 1);
	BOOST_CHECK_CLOSE(mortar(1, 1), 0.75, 1);
	BOOST_CHECK_SMALL(mortar(1, 2), 1e-7);
	BOOST_CHECK_SMALL(mortar(1, 3), 1e-7);
	BOOST_CHECK_SMALL(mortar(1, 4), 1e-7);
	BOOST_CHECK_SMALL(mortar(1, 5), 1e-7);

	BOOST_CHECK_SMALL(mortar(2, 0), 1e-7);
	BOOST_CHECK_CLOSE(mortar(2, 1), 0.75, 1);
	BOOST_CHECK_CLOSE(mortar(2, 2), 0.25, 1);
	BOOST_CHECK_SMALL(mortar(2, 3), 1e-7);
	BOOST_CHECK_SMALL(mortar(2, 4), 1e-7);
	BOOST_CHECK_SMALL(mortar(2, 5), 1e-7);

	BOOST_CHECK_SMALL(mortar(3, 0), 1e-7);
	BOOST_CHECK_CLOSE(mortar(3, 1), 0.125, 1);
	BOOST_CHECK_CLOSE(mortar(3, 2), 0.7917, 1);
	BOOST_CHECK_CLOSE(mortar(3, 3), 0.0833, 1);
	BOOST_CHECK_SMALL(mortar(3, 4), 1e-7);
	BOOST_CHECK_SMALL(mortar(3, 5), 1e-7);

	BOOST_CHECK_SMALL(mortar(4, 0), 1e-7);
	BOOST_CHECK_SMALL(mortar(4, 1), 1e-7);
	BOOST_CHECK_CLOSE(mortar(4, 2), 0.5, 5);
	BOOST_CHECK_CLOSE(mortar(4, 3), 0.5, 5);
	BOOST_CHECK_SMALL(mortar(4, 4), 1e-7);
	BOOST_CHECK_SMALL(mortar(4, 5), 1e-7);

	BOOST_CHECK_SMALL(mortar(5, 0), 1e-7);
	BOOST_CHECK_SMALL(mortar(5, 1), 1e-7);
	BOOST_CHECK_CLOSE(mortar(5, 2), 0.0833, 5);
	BOOST_CHECK_CLOSE(mortar(5, 3), 0.7917, 1);
	BOOST_CHECK_CLOSE(mortar(5, 4), 0.125, 5);
	BOOST_CHECK_SMALL(mortar(5, 5), 1e-7);

	BOOST_CHECK_SMALL(mortar(6, 0), 1e-7);
	BOOST_CHECK_SMALL(mortar(6, 1), 1e-7);
	BOOST_CHECK_SMALL(mortar(6, 2), 1e-7);
	BOOST_CHECK_CLOSE(mortar(6, 3), 0.25, 5);
	BOOST_CHECK_CLOSE(mortar(6, 4), 0.75, 5);
	BOOST_CHECK_SMALL(mortar(6, 5), 1e-7);

	BOOST_CHECK_SMALL(mortar(7, 0), 1e-7);
	BOOST_CHECK_SMALL(mortar(7, 1), 1e-7);
	BOOST_CHECK_SMALL(mortar(7, 2), 1e-7);
	BOOST_CHECK_SMALL(mortar(7, 3), 1e-7);
	BOOST_CHECK_CLOSE(mortar(7, 4), 0.75, 1);
	BOOST_CHECK_CLOSE(mortar(7, 5), 0.25, 1);

	BOOST_CHECK_SMALL(mortar(8, 0), 1e-7);
	BOOST_CHECK_SMALL(mortar(8, 1), 1e-7);
	BOOST_CHECK_SMALL(mortar(8, 2), 1e-7);
	BOOST_CHECK_SMALL(mortar(8, 3), 1e-7);
	BOOST_CHECK_SMALL(mortar(8, 4), 1e-7);
	BOOST_CHECK_CLOSE(mortar(8, 5), 1, 1);
}

template<typename T>  void checkIdentity(const std::vector<std::pair<int, std::vector<std::pair<int, double>>>>& connections)
{
	for (auto&& conn : connections) {
		auto id = conn.first;
		auto&& vec = conn.second;
		BOOST_CHECK_EQUAL(vec.size(), 1);
		auto&& c = vec[0];
		BOOST_CHECK_CLOSE(c.second, 1.0, 1e-8);
	}
}

void display_connections(const std::vector<std::pair<int, std::vector<std::pair<int, double>>>>& connections)
{
	for (auto&& conn : connections) {
		std::cout << "[" << conn.first << "]" << " -> ";
		for (auto&& pair : conn.second) {
			std::cout << "[" << pair.second << "]*" << pair.first << "+";
		}
		std::cout << "\n";
	}
}

BOOST_AUTO_TEST_SUITE(mortar_connector)

BOOST_AUTO_TEST_CASE(ryplpatzak)
{
  squanch::Model<double> model;
  auto&& rect1 = create_rectangular_area<double>(model, 0, 0, 10, 20, 0, 3, {}, { 0.33, 0.33, 0.66, 0.66 });
  auto&& rect2 = create_rectangular_area<double>(model, 10, 0, 10, 20, 0, 2, {}, { 0.33, 0.66 });

  ozp::quadrature::Gaussian<6> g;

  auto s1 = squanch::NurbsSkeleton<double, 1>(rect1, 1, squanch::LineEnds::End);
  auto s2 = squanch::NurbsSkeleton<double, 1>(rect2, 1, squanch::LineEnds::Start);

	std::vector<std::pair<int, std::vector<std::pair<int, double>>>> connections;

  squanch::connect_with_mortar(model, s2, s1, g, connections);

	//checkRpMortar(mortar);
}

BOOST_AUTO_TEST_CASE(ryplpatzakReverse)
{
	squanch::Model<double> model;
	auto&& rect1 = create_rectangular_area<double>(model, 0, 0, 10, 20, 0, 3, {}, { 0.33, 0.33, 0.66, 0.66 });
	auto&& rect2 = create_rectangular_area<double>(model, 10, 0, 10, 20, 0, 2, {}, { 0.33, 0.66,});

	ozp::quadrature::Gaussian<6> g;

	auto s1 = squanch::NurbsSkeleton<double, 1>(rect1, 1, squanch::LineEnds::End);
	auto s2 = squanch::NurbsSkeleton<double, 1>(rect2, 1, squanch::LineEnds::Start);

	std::vector<std::pair<int, std::vector<std::pair<int, double>>>> connections;

	squanch::connect_with_mortar(model, s1, s2, g, connections);

}

BOOST_AUTO_TEST_CASE(dornish1_2d)
{
	squanch::Model<double> model;
	auto&& rect1 = create_rectangular_area<double>(model, 0, 0, 5, 1, 0, 2, {}, {});
	auto&& rect2 = create_rectangular_area<double>(model, 5, 0, 5, 1, 0, 2, {}, {});

	auto set_point = [&model](squanch::Nurbs<double, 2>& rect, int cpix, int cpiy, double x, double y, double w)
	{
		auto&& p = model.get_point(rect.cpi()(cpix, cpiy));
		p[0] = x / w;
		p[1] = y / w;
		p[2] = 0;
		p[3] = w;
	};

	set_point(rect1, 1, 0, 5, 0, 1);
	set_point(rect1, 1, 1, 5.5, 0.3, 1.2);
	set_point(rect1, 1, 2, 5.7, 0.6, 1.5);
	set_point(rect1, 1, 3, 5, 1, 1);

	set_point(rect2, 0, 0, 5, 0, 1);
	set_point(rect2, 0, 1, 5.5, 0.3, 1.2);
	set_point(rect2, 0, 2, 5.7, 0.6, 1.5);
	set_point(rect2, 0, 3, 5, 1, 1);

	h_refine_elements(model, rect1, 1, 1);
	h_refine_elements(model, rect2, 1, 1);

	ozp::quadrature::Gaussian<4> g;
	std::vector<std::pair<int, std::vector<std::pair<int, double>>>> connections1, connections2;

	auto s1 = squanch::NurbsSkeleton<double, 1>(rect1, 1, squanch::LineEnds::End);
	auto s2 = squanch::NurbsSkeleton<double, 1>(rect2, 1, squanch::LineEnds::Start);

	squanch::connect_with_mortar_two_way(model, s1, s2, g, connections1, connections2);

	checkIdentity<double>(connections1);
	checkIdentity<double>(connections2);
}


BOOST_AUTO_TEST_CASE(dornish1_3d)
{
	squanch::Model<double> model;
	auto&& rect1 = create_rectangular_area<double>(model, 0, 0, 5, 1, 0, 2, {}, {});
	auto&& rect2 = create_rectangular_area<double>(model, 5, 0, 5, 1, 0, 2, {}, {});

	auto set_point = [&model](squanch::Nurbs<double, 2>& rect, int cpix, int cpiy, double x, double y, double w)
	{
		auto&& p = model.get_point(rect.cpi()(cpix, cpiy));
		p[0] = x / w;
		p[1] = y / w;
		p[2] = 0;
		p[3] = w;
	};

	set_point(rect1, 1, 0, 5, 0, 1);
	set_point(rect1, 1, 1, 5.5, 0.3, 1.2);
	set_point(rect1, 1, 2, 5.7, 0.6, 1.5);
	set_point(rect1, 1, 3, 5, 1, 1);

	set_point(rect2, 0, 0, 5, 0, 1);
	set_point(rect2, 0, 1, 5.5, 0.3, 1.2);
	set_point(rect2, 0, 2, 5.7, 0.6, 1.5);
	set_point(rect2, 0, 3, 5, 1, 1);

	Eigen::Vector4d extrude_vec(0, 0, 1, 0);
	auto&& prism1 = extrude(model, rect1, extrude_vec);
	auto&& prism2 = extrude(model, rect2, extrude_vec);

	squanch::p_refine(model, prism1, { { 2, 0, 2 } });
	squanch::p_refine(model, prism2, { { 2, 0, 2 } });

	h_refine_elements(model, prism1, 10, 8, 8);
	h_refine_elements(model, prism2, 10, 8, 8);

	ozp::quadrature::Gaussian<4> g;
	std::vector<std::pair<int, std::vector<std::pair<int, double>>>> connections1, connections2;

	auto s1 = squanch::NurbsSkeleton<double, 2>(prism1, { 1, 2 }, squanch::LineEnds::End);
	auto s2 = squanch::NurbsSkeleton<double, 2>(prism2, { 1, 2 }, squanch::LineEnds::Start);
	
	int n0 = prism1.curves()[0].n();
	int n1 = prism1.curves()[1].n();
	int n2 = prism1.curves()[2].n();

	squanch::connect_with_mortar_two_way(model, s1, s2, g, connections1, connections2);

	checkIdentity<double>(connections1);
	checkIdentity<double>(connections2);
}

BOOST_AUTO_TEST_CASE(dornish1_3d_nonconfirming)
{
	squanch::Model<double> model;
	auto&& rect1 = create_rectangular_area<double>(model, 0, 0, 5, 1, 0, 2, {}, {});
	auto&& rect2 = create_rectangular_area<double>(model, 5, 0, 5, 1, 0, 2, {}, {});

	auto set_point = [&model](squanch::Nurbs<double, 2>& rect, int cpix, int cpiy, double x, double y, double w)
	{
		auto&& p = model.get_point(rect.cpi()(cpix, cpiy));
		p[0] = x / w;
		p[1] = y / w;
		p[2] = 0;
		p[3] = w;
	};

	set_point(rect1, 1, 0, 5, 0, 1);
	set_point(rect1, 1, 1, 5.5, 0.3, 1.2);
	set_point(rect1, 1, 2, 5.7, 0.6, 1.5);
	set_point(rect1, 1, 3, 5, 1, 1);

	set_point(rect2, 0, 0, 5, 0, 1);
	set_point(rect2, 0, 1, 5.5, 0.3, 1.2);
	set_point(rect2, 0, 2, 5.7, 0.6, 1.5);
	set_point(rect2, 0, 3, 5, 1, 1);

	Eigen::Vector4d extrude_vec(0, 0, 1, 0);
	auto&& prism1 = extrude(model, rect1, extrude_vec);
	auto&& prism2 = extrude(model, rect2, extrude_vec);

	squanch::p_refine(model, prism1, { { 0, 0, 2 } });
	squanch::p_refine(model, prism2, { { 0, 0, 2 } });

	h_refine_elements(model, prism1, 10, 32, 32);
	h_refine_elements(model, prism2, 10, 48, 48);

	ozp::quadrature::Gaussian<4> g;
	std::vector<std::pair<int, std::vector<std::pair<int, double>>>> connections1, connections2;

	auto s1 = squanch::NurbsSkeleton<double, 2>(prism1, { 1, 2 }, squanch::LineEnds::End);
	auto s2 = squanch::NurbsSkeleton<double, 2>(prism2, { 1, 2 }, squanch::LineEnds::Start);

	int n0 = prism1.curves()[0].n();
	int n1 = prism1.curves()[1].n();
	int n2 = prism1.curves()[2].n();

	squanch::connect_with_mortar_two_way(model, s2, s1, g, connections1, connections2);
}

BOOST_AUTO_TEST_CASE(same_rect)
{
  squanch::Model<double> model;
  auto&& rect1 = create_rectangular_area<double>(model, 0, 0, 10, 20, 0, 3, {}, { 0.33, 0.33, 0.66, 0.66 });
  auto&& rect2 = create_rectangular_area<double>(model, 0, 0, 10, 20, 0, 3, {}, { 0.33, 0.33, 0.66, 0.66 });

  ozp::quadrature::Gaussian<6> g;

	std::vector<std::pair<int, std::vector<std::pair<int, double>>>> connections1, connections2;
  squanch::connect_with_mortar_two_way(model, rect1, rect2, g, connections1, connections2);

	checkIdentity<double>(connections1);
	checkIdentity<double>(connections2);
}

BOOST_AUTO_TEST_SUITE_END()
