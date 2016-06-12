#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "Mesh.h"
#define WIDTH 64
#define HEIGHT 64

int gridX = 600;
int gridY = 600;
int gridZ = 600;

const double fovy = 50.;
const double clipNear = .01;
const double clipFar = 1000.;
double x = 0;
double y = 0;
double z = 3.5;
int currTri = 0;
bool texture = true;

std::string path = "/Users/rohansawhney/Desktop/developer/C++/conformal-parameterization/bunnyhead.obj";

Mesh mesh;
bool success = true;
GLubyte checkerboard[WIDTH][HEIGHT][4];
GLuint texName;

void printInstructions()
{
    std::cerr << "' ': toggle texture/wireframe\n"
              << "↑/↓: move in/out\n"
              << "→/←: change triangle in wireframe mode\n"
              << "w/s: move up/down\n"
              << "a/d: move left/right\n"
              << "escape: exit program\n"
              << std::endl;
}

void buildCheckerboard()
{
    for (int i = 0; i < WIDTH; i++) {
        for (int j = 0; j < HEIGHT; j++) {
            int c = (((i & 0x8) == 0) ^ ((j & 0x8) == 0)) * 255;
            
            if (c == 255) checkerboard[i][j][0] = 0;
            else checkerboard[i][j][0] = (GLubyte) c;
            
            if (c == 255) checkerboard[i][j][1] = 0;
            else checkerboard[i][j][1] = (GLubyte) c;
            
            checkerboard[i][j][2] = (GLubyte) c;
            checkerboard[i][j][3] = (GLubyte) 255;
        }
    }
}

// Initialize OpenGL
void init()
{
    glClearColor (0.1, 0.1, 0.1, 0.0);
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    
    buildCheckerboard();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    
    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                    GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                    GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH,
                 HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                 checkerboard);
}

void drawTexture()
{
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBindTexture(GL_TEXTURE_2D, texName);
    
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        
        if (f->isBoundary()) continue;
        
        // draw flattened
        glBegin(GL_TRIANGLES);
        HalfEdgeCIter he = f->he;
        do {
            glTexCoord2d(he->vertex->uv.x(), he->vertex->uv.y());
            glVertex3d(he->vertex->uv.x() - 0.55,
                       he->vertex->uv.y(),
                       0);
            
            he = he->next;
        } while (he != f->he);
        glEnd();
        
        // draw positions
        glBegin(GL_TRIANGLES);
        do {
            glTexCoord2d(he->vertex->uv.x(), he->vertex->uv.y());
            glVertex3d(he->vertex->position.x() + 1.0,
                       he->vertex->position.y(),
                       he->vertex->position.z());
            
            he = he->next;
        } while (he != f->he);
        glEnd();
    }
}

bool faceContainsEdge(const int vId1, const int vId2)
{
    int c = 0;
    HalfEdgeCIter he = mesh.faces[currTri].he;
    do {
        if (he->vertex->index == vId1 || he->vertex->index == vId2) c++;
        
        he = he->next;
        
    } while (he != mesh.faces[currTri].he);
    
    if (c == 2) return true;
    return false;
}

void drawWireframe()
{
    glDisable(GL_TEXTURE_2D);
    
    glBegin(GL_LINES);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        
        VertexCIter v1 = e->he->vertex;
        VertexCIter v2 = e->he->flip->vertex;
        
        if (faceContainsEdge(v1->index, v2->index)) glColor4f(0.0, 1.0, 0.0, 0.6);
        else glColor4f(0.0, 0.0, 1.0, 0.6);
        
        glVertex3d(v1->uv.x() - 0.55,
                   v1->uv.y(),
                   0);
        glVertex3d(v2->uv.x() - 0.55,
                   v2->uv.y(),
                   0);
            
        glVertex3d(v1->position.x() + 1.0,
                   v1->position.y(),
                   v1->position.z());
        glVertex3d(v2->position.x() + 1.0,
                   v2->position.y(),
                   v2->position.z());
    }
    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = (double)viewport[2] / (double)viewport[3];
    gluPerspective(fovy, aspect, clipNear, clipFar);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 1, z, x, y, 0, 0, 1, 0);
    
    if (success) {
        if (texture) drawTexture();
        else drawWireframe();
    }

    glutSwapBuffers();
}

void keyboard(unsigned char key, int x0, int y0)
{
    switch (key) {
        case 27 :
            exit(0);
        case ' ':
            texture = !texture;
            break;
        case 'a':
            x -= 0.03;
            break;
        case 'd':
            x += 0.03;
            break;
        case 'w':
            y += 0.03;
            break;
        case 's':
            y -= 0.03;
            break;
    }
    
    glutPostRedisplay();
}

void special(int i, int x0, int y0)
{
    switch (i) {
        case GLUT_KEY_UP:
            z -= 0.03;
            break;
        case GLUT_KEY_DOWN:
            z += 0.03;
            break;
        case GLUT_KEY_LEFT:
            if (!texture) {
                currTri --;
                if (currTri < 0) currTri = (int)mesh.faces.size()-1;
            }
            break;
        case GLUT_KEY_RIGHT:
            if (!texture) {
                currTri ++;
                if (currTri == (int)mesh.faces.size()) currTri = 0;
            }
            break;
    }
    
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    success = mesh.read(path);
    mesh.parameterize();
    
    printInstructions();
    glutInitWindowSize(gridX, gridY);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    glutCreateWindow("Conformal Parameterization");
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMainLoop();
    
    return 0;
}
