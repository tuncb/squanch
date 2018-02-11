#pragma once

template <typename T> squanch::Nurbs<T, 2>& create_rectangular_area(squanch::Model<double>& model, T bx, T by, T dx, T dy, unsigned int px, unsigned py, const std::vector<T>& hx, const std::vector<T>& hy)
{
  Eigen::Vector3d bottom_left(bx, by, 0);
  Eigen::Vector3d bottom_right(bx + dx, by, 0);

  squanch::Line<double, Eigen::Vector3d> line(bottom_left, bottom_right);
  auto&& bottom = line.to_nurbs(model);
  squanch::PointContainerTypes<T>::point_type dir;
  dir << 0, dy, 0, 1;
  auto&& nurbs = squanch::extrude(model, bottom, dir);
  std::array<unsigned int, 2> p_ref = { px, py };
  squanch::p_refine(model, nurbs, p_ref);
  squanch::h_refine(model, nurbs, hx, hy);

  return nurbs;
}