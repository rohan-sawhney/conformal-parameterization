#ifndef RENDER_DATA_H
#define RENDER_DATA_H

#include "Shader.h"
#include "Mesh.h"

struct GLVertex {
    Eigen::Vector3f position;
    Eigen::Vector3f normal;
    Eigen::Vector3f barycenter;
    Eigen::Vector3f color;
    Eigen::Vector2f uv;
};

struct GLPickVertex {
    Eigen::Vector3f position;
    Eigen::Vector3f color;
};

class GLMesh {
public:
    // constructor
    GLMesh(Mesh& mesh0);
    
    // destructor
    ~GLMesh();
    
    // setup
    void setup(const std::vector<Eigen::Vector3f>& colors);
    
    // update
    void update(const std::vector<Eigen::Vector3f>& colors);
    
    // draw
    void draw(Shader& shader) const;
    
    // draw
    void drawPick(Shader& shader) const;
    
    // member variables
    std::vector<GLVertex> vertices;
    std::vector<GLPickVertex> pickVertices;
    Mesh& mesh;

private:
    // fills buffers
    void fillBuffers(const std::vector<Eigen::Vector3f>& colors);
    
    // fills pick buffers
    void fillPickBuffers();
    
    // member variables
    GLuint vao, pickVao;
    GLuint vbo, pickVbo;
};

#endif
