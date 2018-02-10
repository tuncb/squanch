#pragma once
#include <iostream>
#include <squanch/Nurbs.h>
#include <squanch/PointContainerTypes.h>

template <typename T> void print_cps(const squanch::Model<T>& model, const squanch::Nurbs<T, 1>& nrb)
{
  auto arr = nrb.cpi().shape();
  for (int i = 0; i < arr[0]; ++i) {
    auto id = nrb.cpi()(i);
    std::cout << "[" << i << "," << "][" << id << "] = ";
    auto& p = model.get_point(id);
    std::cout << p.transpose() << std::endl;
  }
}

template <typename T> void print_cps(const squanch::Model<T>& model, const squanch::Nurbs<T, 2>& nrb)
{
  auto arr = nrb.cpi().shape();
  for (int i = 0; i < arr[1]; ++i) {
    for (int j = 0; j < arr[0]; ++j) {
      auto id = nrb.cpi()(j, i);
      std::cout << "[" << j << "," << i << "][" << id << "] = ";
			auto& p = model.get_point(id);
      std::cout << p.transpose() << std::endl;
    }
  }
}

template <typename T> void print_cps(const squanch::Model<T>& model, const squanch::Nurbs<T, 3>& nrb)
{
  auto arr = nrb.cpi().shape();
  for (int i = 0; i < arr[2]; ++i) {
    for (int j = 0; j < arr[1]; ++j) {
      for (int k = 0; k < arr[0]; ++k) {
        auto id = nrb.cpi()(k, j, i);
        std::cout << "[" << k << "," << j << "," << i << "][" << id << "] = ";
				auto& p = model.get_point(id);
        std::cout << p.transpose() << std::endl;
      }
    }
  }
}