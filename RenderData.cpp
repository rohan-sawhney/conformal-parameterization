#include "RenderData.h"

GLMesh::GLMesh(Mesh& mesh0):
mesh(mesh0),
vao(0),
vbo(0)
{
    
}

GLMesh::~GLMesh()
{
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void GLMesh::fillBuffers(const std::vector<Eigen::Vector3f>& colors)
{
    vertices.resize(3*(mesh.faces.size()-1));
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int i = 0;
            int fIdx = 3*f->index;
            
            HalfEdgeCIter h = f->he;
            do {
                VertexCIter v = h->vertex;
                
                // set vertex
                vertices[fIdx+i].position = h->vertex->position.cast<float>();
                vertices[fIdx+i].normal = v->normal().cast<float>();
                vertices[fIdx+i].color = colors[f->index];
                vertices[fIdx+i].uv = v->uv.cast<float>();
                i++;
                
                h = h->next;
            } while (h != f->he);
            
            // set barycenters
            vertices[fIdx+0].barycenter = Eigen::Vector3f(1.0, 0.0, 0.0);
            vertices[fIdx+1].barycenter = Eigen::Vector3f(0.0, 1.0, 0.0);
            vertices[fIdx+2].barycenter = Eigen::Vector3f(0.0, 0.0, 1.0);
        }
    }
}

Eigen::Vector3f elementColor(unsigned int idx)
{
    return Eigen::Vector3f(((idx & 0x000000FF) >>  0)/255.0,
                           ((idx & 0x0000FF00) >>  8)/255.0,
                           ((idx & 0x00FF0000) >> 16)/255.0);
}

void GLMesh::fillPickBuffers()
{
    int tris = (int)mesh.faces.size()-1;
    
    // add faces
    int i = 0;
    pickVertices.resize(3*tris);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (!f->isBoundary()) {
            int j = 0;
            Eigen::Vector3f color = elementColor(f->index);
            
            HalfEdgeCIter h = f->he;
            do {
                // set vertex
                pickVertices[3*i+j].position = h->vertex->position.cast<float>();
                pickVertices[3*i+j].color = color;
                j++;
                
                h = h->next;
            } while (h != f->he);
            i++;
        }
    }
    
    // add vertex faces
    int verts = 0;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        verts += v->degree();
        if (v->isBoundary()) verts -= 1;
    }
    
    pickVertices.resize(3*(tris + verts));
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        Eigen::Vector3d n = 0.0015*v->normal();
        Eigen::Vector3f p1 = (v->position + n).cast<float>();
        Eigen::Vector3f color = elementColor(tris + v->index);
        
        HalfEdgeCIter h = v->he;
        do {
            if (!h->onBoundary) {
                Eigen::Vector3f p2 = (h->next->vertex->position + n).cast<float>();
                Eigen::Vector3f p3 = (h->next->next->vertex->position + n).cast<float>();
                
                pickVertices[3*i+0].position = p1;
                pickVertices[3*i+1].position = p1 + 0.125*(p2 - p1);
                pickVertices[3*i+2].position = p1 + 0.125*(p3 - p1);
                pickVertices[3*i+0].color = pickVertices[3*i+1].color = pickVertices[3*i+2].color = color;
                i++;
            }
            h = h->flip->next;
        } while (h != v->he);
    }
}

void GLMesh::setup(const std::vector<Eigen::Vector3f>& colors)
{
    fillBuffers(colors);
    fillPickBuffers();
    
    // bind vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // bind vbo
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizei)vertices.size() * sizeof(GLVertex), &vertices[0], GL_STATIC_DRAW);

    // set vertex attribute pointers for vbo data
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)0);
    
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, normal));
    
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, barycenter));
    
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, color));
    
    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, sizeof(GLVertex), (GLvoid*)offsetof(GLVertex, uv));
    
    // bind pick vao
    glGenVertexArrays(1, &pickVao);
    glBindVertexArray(pickVao);
    
    // bind pick vbo
    glGenBuffers(1, &pickVbo);
    glBindBuffer(GL_ARRAY_BUFFER, pickVbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizei)pickVertices.size() * sizeof(GLPickVertex),
                 &pickVertices[0], GL_STATIC_DRAW);
    
    // set vertex attribute pointers for pick vbo data
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLPickVertex), (GLvoid*)0);
    
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GLPickVertex), (GLvoid*)offsetof(GLPickVertex, color));
    
    // unbind vao
    glBindVertexArray(0);
}

void GLMesh::update(const std::vector<Eigen::Vector3f>& colors)
{
    fillBuffers(colors);
    fillPickBuffers();
    
    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizei)vertices.size() * sizeof(GLVertex), &vertices[0]);
    
    // bind pick vbo
    glBindBuffer(GL_ARRAY_BUFFER, pickVbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizei)pickVertices.size() * sizeof(GLPickVertex), &pickVertices[0]);
}

void GLMesh::draw(Shader& shader) const
{
    shader.use();    
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)vertices.size());
    glBindVertexArray(0);
}

void GLMesh::drawPick(Shader& shader) const
{
    shader.use();
    glBindVertexArray(pickVao);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)pickVertices.size());
    glBindVertexArray(0);
}
