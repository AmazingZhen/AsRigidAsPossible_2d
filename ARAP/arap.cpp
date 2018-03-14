#include "arap.h"

#include <iostream>

void ARAP::registration(const std::vector<float>& vertices, const std::vector<unsigned int> &indices)
{
	int row = indices.size(), col = vertices.size() / 3;

	// Min ||(v1' - v0') - (v1 - v0)||
	A_error = Eigen::MatrixXf::Zero(row, col);
	b_x_error = Eigen::MatrixXf::Zero(row, 1);
	b_y_error = Eigen::MatrixXf::Zero(row, 1);

	std::vector<std::vector<bool>> edgeUsed(vertices.size(), std::vector<bool>(vertices.size(), false));
	int curRow = 0;
	for (int i = 0; i < indices.size(); i += 3) {
		for (int j = i; j < i + 3; j++) {
			int vStart = indices[j], vEnd = indices[(j + 1 == i + 3 ? i : j + 1)];

			A_error(curRow, vStart) = -1;
			A_error(curRow, vEnd) = 1;

			b_x_error(curRow, 0) = vertices[vEnd * 3] - vertices[vStart * 3];
			b_y_error(curRow, 0) = vertices[vEnd * 3 + 1] - vertices[vStart * 3 + 1];
			
			curRow++;
		}
	}

	//printf("A_Error:\n");
	//std::cout << A_error << std::endl;
	//std::cout << b_x_error << std::endl;
	//std::cout << b_y_error << std::endl;
	// A * v'x = bx and A * v'y = by
}

// Add constraints for select vetices.
void ARAP::compilation(const std::unordered_set<unsigned int> vSelected)
{
	A_constraint = Eigen::MatrixXf::Zero(vSelected.size(), A_error.cols());
	b_x_constraint = Eigen::MatrixXf::Zero(vSelected.size(), 1);
	b_y_constraint = Eigen::MatrixXf::Zero(vSelected.size(), 1);

	int curRow = 0;
	w = 100.f;  // w increases, constraints increases.
	for (auto iter = vSelected.begin(); iter != vSelected.end(); iter++) {
		int index = *iter;
		A_constraint(curRow, index) = w;
		curRow++;

		//printf("choose index: %d\n", index);
	}

	//printf("A_constraint:\n");
	//std::cout << A_constraint << std::endl;
}

void ARAP::solve(std::vector<float>& vertices, const std::unordered_set<unsigned int> vSelected)
{
	// Update constraint vertices
	int curRow = 0;
	for (auto iter = vSelected.begin(); iter != vSelected.end(); iter++) {
		int index = *iter;
		b_x_constraint(curRow, 0) = w * vertices[index * 3];
		b_y_constraint(curRow, 0) = w * vertices[index * 3 + 1];
		curRow++;
	}

	Eigen::MatrixXf A(A_error.rows() + A_constraint.rows(), A_error.cols());
	A << A_error, A_constraint;

	//printf("A:\n");
	//std::cout << A << std::endl;

	Eigen::MatrixXf b_x(b_x_error.rows() + b_x_constraint.rows(), b_x_error.cols());
	b_x << b_x_error, b_x_constraint;
	
	Eigen::MatrixXf b_y(b_y_error.rows() + b_y_constraint.rows(), b_y_error.cols());
	b_y << b_y_error, b_y_constraint;

	// v = (AT * A)^-1 * AT * b
	Eigen::MatrixXf pinv = (A.transpose() * A).completeOrthogonalDecomposition().pseudoInverse();

	Eigen::MatrixXf v_x = pinv * A.transpose() * b_x;
	Eigen::MatrixXf v_y = pinv * A.transpose() * b_y;

	for (int i = 0; i < v_x.rows(); i++) {
		//printf("index %d: (%f, %f), (%f, %f)\n", i, vertices[i * 3], vertices[i * 3 + 1], v_x(i, 0), v_y(i, 0));
		vertices[i * 3] = v_x(i, 0);
		vertices[i * 3 + 1] = v_y(i, 0);
	}

	return;
}
