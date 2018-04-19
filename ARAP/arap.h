#pragma once

#include <vector>
#include <set>
#include <map>

#include <Eigen/Dense>

enum Solver_mode { base, similarity, similarity_scale };

class ARAP {
public:
	ARAP();

	void registration(const std::vector<float> &vertices,
		const std::vector<unsigned int> &indices);

	void compilation(const std::set<unsigned int> vSelected);

	void solve(std::vector<float> &vertices, const std::vector<unsigned int>& indices, const std::set<unsigned int> vSelected);

private:
	void registration_base(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void registration_similarity(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void registration_similarity_scale(const std::vector<float>& vertices, const std::vector<unsigned int>& indices);

	void solve_base(std::vector<float> &vertices, const std::set<unsigned int> vSelected);

	void solve_similarity(std::vector<float> &vertices, const std::set<unsigned int> vSelected);

	void solve_similarity_scale(std::vector<float> &vertices,
		const std::vector<unsigned int>& indices,
		const std::set<unsigned int> vSelected);

	Solver_mode mode;

	Eigen::MatrixXf A1_error;
	Eigen::MatrixXf b1_error;

	Eigen::MatrixXf A1_constraint;
	Eigen::MatrixXf b1_constraint;

	std::vector<Eigen::MatrixXf> G_all;

	Eigen::MatrixXf A2_error;
	Eigen::MatrixXf b2_error;

	Eigen::MatrixXf A1;
	Eigen::MatrixXf A2;

	Eigen::MatrixXf b1;
	Eigen::MatrixXf b2;
	Eigen::MatrixXf b2_error_transformed;

	Eigen::LLT<Eigen::MatrixXf> A1T_A1_LLT;
	Eigen::LLT<Eigen::MatrixXf> A2T_A2_LLT;

	Eigen::MatrixXf A1_errorT_b1_error;

	float w;

	std::vector<std::vector<bool>> edgeUsed;
	std::map<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>> edgeTableForTriangleIndex;
};