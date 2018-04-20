#include "render.h"

GLuint vbo, ebo;
float gWidth = 1000, gHeight = 1000;
float gScale = 1.0f;
float gX, gY;
bool gLeftButtonPressed = false;

const unsigned int nRowLen = 6;

Render *r = 0;

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_W && action == GLFW_PRESS) {
		
	}
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	gScale += yoffset / 10.f ;
}

int selectIndex = -1;

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (action == GLFW_PRESS) {
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			gLeftButtonPressed = true;

			float world_x = (float((gX - gWidth / 2) / gWidth) * 2) / gScale;
			float world_y = (float(-(gY - gHeight / 2) / gHeight) * 2) / gScale;

			selectIndex = r->findClosestVertexIndex(world_x, world_y);
		}
		else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
			float world_x = (float((gX - gWidth / 2) / gWidth) * 2) / gScale;
			float world_y = (float(-(gY - gHeight / 2) / gHeight) * 2) / gScale;

			r->selectVertex(world_x, world_y);
		}
	}
	else {
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			gLeftButtonPressed = false;
		}
	}
}

void cursor_position_callback(GLFWwindow* window, double x, double y)
{
	gX = x;
	gY = y;

	if (gLeftButtonPressed) {
		float world_x = (float((gX - gWidth / 2) / gWidth) * 2) / gScale;
		float world_y = (float(-(gY - gHeight / 2) / gHeight) * 2) / gScale;

		r->moveVertexSelected(world_x, world_y);
	}

	return;
}

bool Render::init()
{
	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		return false;
	}

	glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	//glfwWindowHint(GLFW_CONTEXT_PROFILE, GLFW_OPENGL_CORE_PROFILE); //We don't want the old OpenGL

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(gWidth, gHeight, "ARAP", NULL, NULL);

	if (window == NULL)
	{
		fprintf(stderr, "Failed to open GLFW window\n");
		glfwTerminate();
		return false;
	}

	// Initialize GLEW
	glfwMakeContextCurrent(window);
	glewExperimental = true; // Needed in core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return false;
	}

	glfwSetKeyCallback(window, key_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);

	r = this;
	solver = ARAP();

	return true;
}

void Render::start()
{
	if (!init()) {
		return;
	}

	generateSquareMesh();

	solver.registration(vertices, indices);

	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// run while the window is open
	while (!glfwWindowShouldClose(window)) {
		// process pending events
		glfwPollEvents();

		// draw one frame
		render();
	}

	glfwTerminate();
}

void Render::selectVertex(float x, float y)
{
	int minDistIndex = findClosestVertexIndex(x, y);

	if (minDistIndex != -1) {
		auto pos = vSelected.find(minDistIndex);
		if (pos == vSelected.end()) {
			vSelected.insert(minDistIndex);
		}
		else {
			vSelected.erase(pos);
		}

		solver.compilation(vSelected);
		solver.solve(vertices, indices, vSelected);
	}
}

void Render::moveVertexSelected(float x, float y)
{
	// int minDistIndex = findClosestVertexIndex(x, y);
	// printf("selected : %d\n", minDistIndex);

	if (vSelected.find(selectIndex) != vSelected.end()) {
		vertices[selectIndex * 3] = x;
		vertices[selectIndex * 3 + 1] = y;

		solver.solve(vertices, indices, vSelected);
	}

}

void Render::generateSquareMesh()
{
	vertices.clear();
	indices.clear();

	float fYStep = 2.0f / (float)(nRowLen - 1);
	float fXStep = 2.0f / (float)(nRowLen - 1);

	for (unsigned int yi = 0; yi < nRowLen; ++yi) {
		float fY = -1.0f + (float)yi * fYStep;
		for (unsigned int xi = 0; xi < nRowLen; ++xi) {
			float fX = -1.0f + (float)xi * fXStep;
			float v[3] = { fX, fY, 0.f };

			vertices.push_back(v[0]);
			vertices.push_back(v[1]);
			vertices.push_back(v[2]);
		}
	}

	printf("num of v: %d\n", vertices.size() / 3);

	originalVertices = vertices;

	for (unsigned int yi = 0; yi < nRowLen - 1; ++yi) {
		unsigned int nRow1 = yi * nRowLen;
		unsigned int nRow2 = (yi + 1) * nRowLen;

		for (unsigned int xi = 0; xi < nRowLen - 1; ++xi) {
			unsigned int nTri1[3] = { nRow1 + xi, nRow2 + xi + 1, nRow1 + xi + 1 };
			unsigned int nTri2[3] = { nRow1 + xi, nRow2 + xi, nRow2 + xi + 1 };

			indices.push_back(nTri1[0]);
			indices.push_back(nTri1[1]);
			indices.push_back(nTri1[2]);

			indices.push_back(nTri2[0]);
			indices.push_back(nTri2[1]);
			indices.push_back(nTri2[2]);
		}
	}

	/*
	for (int i = 0; i < vertices.size(); i += 3) {
		printf("(%f, %f, %f)\n", vertices[i], vertices[i + 1], vertices[i + 2]);
	}

	for (int i = 0; i < indices.size(); i += 3) {
		printf("{%d, %d, %d}\n", indices[i], indices[i + 1], indices[i + 2]);
	}
	*/

	/*
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, &ebo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	*/
}

void Render::render()
{	
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);
	
	/*
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.f, 0.f, -1.f);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

	glDrawElements(
		GL_TRIANGLES,
		indices.size(),
		GL_UNSIGNED_INT,
		0
	);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	//*/

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(gScale, gScale, 1.0f);


	// Draw triangles;
	glLineWidth(2.0f);
	glColor3f(1.0f, 1.0f, 1.0f);

	for (unsigned int i = 0; i < indices.size(); i += 3) {
		int index[3] = { indices[i], indices[i + 1], indices[i + 2] };

		glBegin(GL_LINE_LOOP);

		glVertex3f(vertices[index[0] * 3], vertices[index[0] * 3 + 1], vertices[index[0] * 3 + 2]);
		glVertex3f(vertices[index[1] * 3], vertices[index[1] * 3 + 1], vertices[index[1] * 3 + 2]);
		glVertex3f(vertices[index[2] * 3], vertices[index[2] * 3 + 1], vertices[index[2] * 3 + 2]);

		glEnd();
	}

	// Draw select vertices
	glColor3f(1.0f, 0.0f, 0.0f);

	for (auto iter = vSelected.begin(); iter != vSelected.end(); iter++) {
		int index = *iter;
		float x = vertices[index * 3], y = vertices[index * 3 + 1], z = vertices[index * 3 + 2];

		glBegin(GL_QUADS);

		glVertex3f(x - 0.02, y - 0.02, 0.f);
		glVertex3f(x + 0.02, y - 0.02, 0.f);
		glVertex3f(x + 0.02, y + 0.02, 0.f);
		glVertex3f(x - 0.02, y + 0.02, 0.f);

		glEnd();
	}

	glfwSwapBuffers(window);
}

int Render::findClosestVertexIndex(float x, float y)
{
	float minDist = INT_MAX;
	int minDistIndex = -1;

	for (int i = 0; i < vertices.size(); i += 3) {
		float dist = 0.f;
		dist += (vertices[i] - x) * (vertices[i] - x);
		dist += (vertices[i + 1] - y) * (vertices[i + 1] - y);

		if (dist < minDist) {
			minDist = dist;
			minDistIndex = i / 3;
		}
	}

	if (minDist < 0.1) {
		return minDistIndex;
	}

	return -1;
}
