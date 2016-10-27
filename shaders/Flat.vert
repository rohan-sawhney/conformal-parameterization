#version 410 core
layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inColor;

out vec3 position;
out vec3 color;

layout (std140) uniform Transform {
    mat4 projection;
    mat4 view;
    mat4 model;
};

void main()
{
    gl_Position = projection * view * model * vec4(inPosition, 1.0);
    position = vec3(model * vec4(inPosition, 1.0));
    color = inColor;
}
