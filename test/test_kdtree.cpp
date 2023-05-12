#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "kdtree.h"
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

TEST_CASE("Test KDTree") {
  srand(time(NULL));
  int n_dim = 5;
  int n_leaf = 4;
  vector<Point> ps;
  for (int i = 0; i < 100; i++) {
    Point p(5, 0);
    for (int i = 0; i < 5; i++) {
      double x = 0 + (rand() % 1000);
      p[i] = x;
    }
    ps.push_back(p);
  }

  KDTree tree = KDTree(n_dim, n_leaf, ps);

  vector<Point> qps = tree.kneighbors(ps[10], 10);
  cout << qps.size() << endl;
  CHECK(qps.size() <= 10);
}