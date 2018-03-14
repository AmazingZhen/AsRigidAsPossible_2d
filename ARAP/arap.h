#pragma once

#include <vector>
#include <unordered_set>

#include <Eigen/Dense>

class ARAP {

public:
	void registration(const std::vector<float> &vertices,
		const std::vector<unsigned int> &indices);

	void compilation(const std::unordered_set<unsigned int> vSelected);

	void solve(std::vector<float> &vertices,
		const std::unordered_set<unsigned int> vSelected);

private:
	Eigen::MatrixXf G;

	Eigen::MatrixXf A_error;
	Eigen::MatrixXf A_constraint;
	Eigen::MatrixXf b_x_error;
	Eigen::MatrixXf b_x_constraint;
	Eigen::MatrixXf b_y_error;
	Eigen::MatrixXf b_y_constraint;
	float w;
};