#pragma once

#include <vector>
#include <unordered_set>

#include <Eigen/Dense>

class ARAP {

public:
	void registration(const std::vector<float> &vertices,
		const std::vector<unsigned int> &indices);

	void registration_base(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void registration_similarity(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void registration_similarity_scale(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void compilation(const std::unordered_set<unsigned int> vSelected);

	void solve(std::vector<float> &vertices,
		const std::unordered_set<unsigned int> vSelected);

private:
	Eigen::MatrixXf A1_error;
	Eigen::MatrixXf b1_error;

	Eigen::MatrixXf A1_constraint;
	Eigen::MatrixXf b1_constraint;

	std::vector<Eigen::MatrixXf> G_all;

	float w;
};