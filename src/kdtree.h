#ifndef KD_TREE_HPP
#define KD_TREE_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

using namespace std;
typedef vector<double> Point;

void print_point(Point &p) {
  for (int i = 0; i < p.size(); i++) {
    cout << p[i] << " ";
  }
  cout << endl;
}

class KDTree {
private:
  int dim;
  int cur_dim;
  int leaf_size;
  Point median;

  // for nearest neighbor search and used by leaf node.
  Point centroid;
  KDTree *left;
  KDTree *right;
  void split();
  vector<pair<double, KDTree *>> get_search_leaves(const Point &q, const int k);
  void get_search_leaves(vector<pair<double, KDTree *>> &inserts,
                         const Point &q);

public:
  vector<Point> points;

  // The constructor for root node
  KDTree(int n_dim, int n_leaf, vector<Point> ps);
  KDTree(int n_dim, int c_dim, int n_leaf, vector<Point> ps);
  ~KDTree();
  vector<Point> kneighbors(const Point &p, const int k);
};

Point calculate_centroid(const vector<Point> &ps, int n_dim) {
  Point p(n_dim, 0);

  for (auto point : ps) {
    for (int i = 0; i < n_dim; i++) {
      p[i] += point[i];
    }
  }

  for (int i = 0; i < n_dim; i++) {
    p[i] /= ps.size();
  }

  return p;
}

KDTree::KDTree(int n_dim, int n_leaf, vector<Point> ps) {
  dim = n_dim;
  cur_dim = 0;
  leaf_size = n_leaf;
  left = nullptr;
  right = nullptr;
  points = ps;

  if (points.size() > leaf_size) {
    split();
  } else {
    centroid = calculate_centroid(points, dim);
    // print_point(centroid);
  }
}

KDTree::KDTree(int n_dim, int c_dim, int n_leaf, vector<Point> ps) {
  dim = n_dim;
  cur_dim = c_dim;
  leaf_size = n_leaf;
  left = nullptr;
  right = nullptr;
  points = ps;

  if (points.size() > leaf_size) {
    split();
  } else {
    centroid = calculate_centroid(points, dim);
    // print_point(centroid);
  }
}

KDTree::~KDTree() {
  delete left;
  delete right;
}

void KDTree::split() {
  sort(points.begin(), points.end(),
       [&](Point a, Point b) { return (a[cur_dim] < b[cur_dim]); });

  int m = points.size() / 2;
  int next_dim = (cur_dim + 1) % dim;
  median = points[m];
  vector<Point> left_points = vector<Point>(points.begin(), points.begin() + m);
  vector<Point> right_points = vector<Point>(points.begin() + m, points.end());

  left = new KDTree(dim, next_dim, leaf_size, left_points);
  right = new KDTree(dim, next_dim, leaf_size, right_points);
  points.clear();
}

double distance(const Point &a, const Point &b, int dim) {
  double sum = 0.0;
  for (int i = 0; i < dim; i++) {
    sum += pow(a[i] - b[i], 2);
  }
  return sqrt(sum);
}

void KDTree::get_search_leaves(vector<pair<double, KDTree *>> &inserts,
                               const Point &q) {
  bool is_leaf = (left == nullptr) && (right == nullptr);

  if (is_leaf) {
    double d = distance(q, centroid, dim);
    inserts.push_back(pair<double, KDTree *>(d, this));
  } else {
    right->get_search_leaves(inserts, q);
    left->get_search_leaves(inserts, q);
  }
}

vector<pair<double, KDTree *>> KDTree::get_search_leaves(const Point &q,
                                                         const int k) {
  vector<pair<double, KDTree *>> leaves;
  bool is_leaf = (left == nullptr) && (right == nullptr);

  if (is_leaf) {
    double d = distance(q, centroid, dim);
    leaves.push_back(pair<double, KDTree *>(d, this));
  } else {
    right->get_search_leaves(leaves, q);
    left->get_search_leaves(leaves, q);
  }

  sort(leaves.begin(), leaves.end(),
       [](pair<double, KDTree *> &a, pair<double, KDTree *> &b) {
         return a.first < b.first;
       });

  return leaves;
}

void query_list_insert(list<pair<double, Point>> &ql, const Point &q, Point p,
                       int dim) {
  double d = distance(q, p, dim);
  if (ql.size() == 0) {
    ql.push_back(pair<double, Point>(d, p));
    return;
  }

  bool has_insert = false;

  for (auto it = ql.begin(); it != ql.end(); it++) {
    double compare = it->first;
    if (d < compare) {
      ql.insert(it, pair<double, Point>(d, p));
      has_insert = true;
      break;
    }
  }

  if (!has_insert) {
    ql.insert(ql.end(), pair<double, Point>(d, p));
  }
}

vector<Point> KDTree::kneighbors(const Point &q, const int k) {
  list<pair<double, Point>> finding;
  vector<pair<double, KDTree *>> leaves = get_search_leaves(q, k);

  auto end = leaves.end();
  if (leaves.size() > k) {
    end = leaves.begin();
    end += k;
  }

  for (auto it = leaves.begin(); it != end; it++) {
    KDTree *tree = it->second;
    for (auto &point : tree->points) {
      query_list_insert(finding, q, point, dim);
    }
  }

  vector<Point> result;
  for (auto it = finding.begin(); it != finding.end(); it++) {
    result.push_back(it->second);
    if (result.size() == k) {
      break;
    }
  }
  return result;
}

#endif