#include "Mesh.h"
#include "RenderData.h"
#include "Camera.h"

#define ESCAPE 27
#define DIGIT_OFFSET 48

int gridX = 600;
int gridY = 600;

GLuint transformUbo;
GLuint lightUbo;

const std::string path = "/Users/rohansawhney/Desktop/developer/C++/conformal-parameterization/objs/bunnyhead.obj";
const std::string shaderPath = "/Users/rohansawhney/Desktop/developer/C++/conformal-parameterization/shaders/";
Shader meshShader(shaderPath);
Shader normalShader(shaderPath);
Shader wireframeShader(shaderPath);
Shader checkerboardShader(shaderPath);

Camera camera;
float lastTime = 0.0;
float dt = 0.0;
float lastX = 0.0, lastY = 0.0;
bool keys[256];
bool firstMouse = true;

Mesh mesh, parameterMesh;
GLMesh glMesh(mesh), glParameterMesh(parameterMesh);
Eigen::Matrix4f transform;

std::vector<Eigen::Vector3f> defaultColors;
std::vector<Eigen::Vector3f> qcErrors;

const Eigen::Vector3f lightPosition(0.0, 3.0, -3.0);
const Eigen::Vector3f lightColor(1.0, 1.0, 1.0);

bool success = true;
bool showCheckerboard = false;
bool showQcError = false;
bool showWireframe = false;
bool showNormals = false;

void setupShaders()
{
    meshShader.setup("Model.vert", "", "Model.frag");
    normalShader.setup("Normal.vert", "Normal.geom", "Normal.frag");
    wireframeShader.setup("Wireframe.vert", "", "Wireframe.frag");
    checkerboardShader.setup("Model.vert", "", "Checkerboard.frag");
}

void setupUniformBlocks()
{
    // 1) generate transform indices
    GLuint meshShaderIndex = glGetUniformBlockIndex(meshShader.program, "Transform");
    GLuint normalShaderIndex = glGetUniformBlockIndex(normalShader.program, "Transform");
    GLuint wireframeShaderIndex = glGetUniformBlockIndex(wireframeShader.program, "Transform");
    GLuint checkerboardShaderIndex = glGetUniformBlockIndex(checkerboardShader.program, "Transform");
    
    // bind
    glUniformBlockBinding(meshShader.program, meshShaderIndex, 0);
    glUniformBlockBinding(normalShader.program, normalShaderIndex, 0);
    glUniformBlockBinding(wireframeShader.program, wireframeShaderIndex, 0);
    glUniformBlockBinding(checkerboardShader.program, checkerboardShaderIndex, 0);
    
    // add transform data
    glGenBuffers(1, &transformUbo);
    glBindBuffer(GL_UNIFORM_BUFFER, transformUbo);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, transformUbo);
    glBufferData(GL_UNIFORM_BUFFER, 3*sizeof(Eigen::Matrix4f), NULL, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    
    // 2) generate light index
    meshShaderIndex = glGetUniformBlockIndex(meshShader.program, "Light");
    
    // bind
    glUniformBlockBinding(meshShader.program, meshShaderIndex, 1);
    
    // add light data
    glGenBuffers(1, &lightUbo);
    glBindBuffer(GL_UNIFORM_BUFFER, lightUbo);
    glBufferData(GL_UNIFORM_BUFFER, 2*sizeof(Eigen::Vector4f), NULL, GL_STATIC_DRAW); // std140 alignment
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, lightUbo);
    
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(Eigen::Vector4f), lightPosition.data());
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(Eigen::Vector4f), sizeof(Eigen::Vector4f), lightColor.data());
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void printInstructions()
{
    std::cerr << "1: scp\n"
              << "2: lscm\n"
              << "3: circle patterns\n"
              << "4: cetm\n"
              << "5: toggle checkerboard\n"
              << "6: toggle quasi conformal error\n"
              << "7: toggle normals\n"
              << "8: toggle wireframe\n"
              << "w/s: move in/out\n"
              << "a/d: move left/right\n"
              << "e/q: move up/down\n"
              << "→/←: update face correspondence\n"
              << "escape: exit program\n"
              << std::endl;
}

void setDefaultColors()
{
    const Eigen::Vector3f blue(0.0, 0.0, 1.0);
    defaultColors.resize(mesh.faces.size());
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        defaultColors[f->index] = blue;
    }
}

void updateErrors()
{
    qcErrors.resize(mesh.faces.size());
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        qcErrors[f->index] = Eigen::Vector3f(f->qcError.x(), f->qcError.y(), f->qcError.z());
    }
}

void init()
{
    // enable depth testing and multisampling
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    
    glClearColor(0.1, 0.1, 0.1, 1.0);
    glClearDepth(1.0);
    glDepthFunc(GL_LEQUAL);
    
    // setup shaders and blocks
    setupShaders();
    setupUniformBlocks();
    
    // read mesh
    success = mesh.read(path) && parameterMesh.read(path);;
    if (success) {
        setDefaultColors();
        glMesh.setup(defaultColors);
        glParameterMesh.setup(defaultColors);
    }
    
    // print instructions
    printInstructions();
}

void uploadCameraTransforms()
{
    // set camera transformations
    glBindBuffer(GL_UNIFORM_BUFFER, transformUbo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4),
                    glm::value_ptr(camera.projectionMatrix(gridX, gridY)));
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4),
                    glm::value_ptr(camera.viewMatrix()));
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    
    // set view position
    meshShader.use();
    glUniform3f(glGetUniformLocation(meshShader.program, "viewPosition"),
                camera.pos.x(), camera.pos.y(), camera.pos.z());
}

void draw(GLMesh& m)
{
    if (showCheckerboard) m.draw(checkerboardShader);
    else m.draw(meshShader);
    if (showNormals) m.draw(normalShader);
    if (showWireframe) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        m.draw(wireframeShader);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
    }
}

void display()
{
    float elapsedTime = glutGet(GLUT_ELAPSED_TIME);
    dt = (elapsedTime - lastTime) / 1000.0;
    lastTime = elapsedTime;
    
    if (success) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // upload camera transforms
        uploadCameraTransforms();
        
        // set model transformation
        glBindBuffer(GL_UNIFORM_BUFFER, transformUbo);
        transform.setIdentity(); transform(0, 3) = -1.0;
        glBufferSubData(GL_UNIFORM_BUFFER, 2*sizeof(Eigen::Matrix4f),
                        sizeof(Eigen::Matrix4f), transform.data());
        
        // draw
        draw(glMesh);
        
        transform.setIdentity(); transform(0, 3) = 1.0;
        glBufferSubData(GL_UNIFORM_BUFFER, 2*sizeof(Eigen::Matrix4f),
                        sizeof(Eigen::Matrix4f), transform.data());
        
        // draw
        draw(glParameterMesh);
        
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        // swap
        glutSwapBuffers();
    }
}

void idle()
{
    glutPostRedisplay();
}

void reset()
{
    glDeleteBuffers(1, &transformUbo);
    glDeleteBuffers(1, &lightUbo);
}

void updateScene(const std::string& algorithm)
{
    // copy uvs
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        parameterMesh.vertices[v->index].position = Eigen::Vector3d(v->uv.x(), v->uv.y(), 0);
        parameterMesh.vertices[v->index].uv = v->uv;
    }
    
    // update errors
    updateErrors();
    
    // update gl meshes
    glMesh.update(showQcError ? qcErrors : defaultColors);
    glParameterMesh.update(showQcError ? qcErrors : defaultColors);
    
    // update title
    std::string title = "Conformal Parameterization: " + algorithm;
    glutSetWindowTitle(title.c_str());
}

void keyboardPressed(unsigned char key, int x, int y)
{
    keys[key] = true;
    
    if (keys[ESCAPE]) {
        reset();
        exit(0);
    
    } else if (keys[DIGIT_OFFSET + 1]) {
        mesh.parameterize(SCP);
        updateScene("SCP");
        
    } else if (keys[DIGIT_OFFSET + 2]) {
        mesh.parameterize(LSCM);
        updateScene("LSCM");
    
    } else if (keys[DIGIT_OFFSET + 3]) {
        mesh.parameterize(CIRCLE_PATTERNS);
        updateScene("CIRCLE PATTERNS");
        
    } else if (keys[DIGIT_OFFSET + 4]) {
        mesh.parameterize(CETM);
        updateScene("CETM");
        
    } else if (keys[DIGIT_OFFSET + 5]) {
        showCheckerboard = !showCheckerboard;
        if (showCheckerboard) showQcError = false;
        glMesh.update(showQcError ? qcErrors : defaultColors);
        glParameterMesh.update(showQcError ? qcErrors : defaultColors);
        
    } else if (keys[DIGIT_OFFSET + 6]) {
        showQcError = !showQcError;
        if (showQcError) showCheckerboard = false;
        glMesh.update(showQcError ? qcErrors : defaultColors);
        glParameterMesh.update(showQcError ? qcErrors : defaultColors);
    
    } else if (keys[DIGIT_OFFSET + 7]) {
        showNormals = !showNormals;
        
    } else if (keys[DIGIT_OFFSET + 8]) {
        showWireframe = !showWireframe;
    
    } else if (keys['a']) {
        camera.processKeyboard(LEFT, dt);
        
    } else if (keys['d']) {
        camera.processKeyboard(RIGHT, dt);
        
    } else if (keys['w']) {
        camera.processKeyboard(FORWARD, dt);
        
    } else if (keys['s']) {
        camera.processKeyboard(BACKWARD, dt);
        
    } else if (keys['e']) {
        camera.processKeyboard(UP, dt);
        
    } else if (keys['q']) {
        camera.processKeyboard(DOWN, dt);
    }
}

void keyboardReleased(unsigned char key, int x, int y)
{
    if (key != ESCAPE) keys[key] = false;
}

void special(int i, int x, int y)
{
    switch (i) {
        case GLUT_KEY_LEFT:
            // TODO
            break;
        case GLUT_KEY_RIGHT:
            // TODO
            break;
    }
}

void mouse(int x, int y)
{
    if (firstMouse) {
        lastX = x;
        lastY = y;
        firstMouse = false;
    }
    
    float dx = x - lastX;
    float dy = lastY - y;
    
    lastX = x;
    lastY = y;
    
    camera.processMouse(dx, dy);
}

int main(int argc, char** argv)
{
    // TODOs:
    // 1) Correspondence
    // 2) Linking
    // 3) Circle Patterns
    // 4) Cetm
    // 5) LBFGS
    // 6) Subdivision
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE | GLUT_MULTISAMPLE);
    glutInitWindowSize(gridX, gridY);
    glutCreateWindow("Conformal Parameterization: SCP");
    
    init();
    mesh.parameterize(SCP);
    updateScene("SCP");
    
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboardPressed);
    glutKeyboardUpFunc(keyboardReleased);
    glutSpecialFunc(special);
    glutMotionFunc(mouse);
    glutMainLoop();
    
    return 0;
}
