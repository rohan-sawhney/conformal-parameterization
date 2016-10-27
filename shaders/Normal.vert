#version 410 core
layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inNormal;

out vec3 normal;

layout (std140) uniform Transform {
    mat4 projection;
    mat4 view;
    mat4 model;
};

void main()
{
    gl_Position = projection * view * model * vec4(inPosition, 1.0);
    normal = normalize(vec3(projection * vec4(mat3(transpose(inverse(view * model))) * inNormal, 1.0)));
}
