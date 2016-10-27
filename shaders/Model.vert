#version 410 core
layout (location = 0) in vec3 inPosition;
layout (location = 1) in vec3 inNormal;
layout (location = 3) in vec3 inColor;
layout (location = 4) in vec2 inTexCoords;

out vec3 position;
out vec3 normal;
out vec3 color;
out vec2 texCoords;

layout (std140) uniform Transform {
    mat4 projection;
    mat4 view;
    mat4 model;
};

void main()
{
    gl_Position = projection * view * model * vec4(inPosition, 1.0);
    position = vec3(model * vec4(inPosition, 1.0));
    normal = normalize(vec3(projection * vec4(mat3(transpose(inverse(view * model))) * inNormal, 1.0)));
    color = inColor;
    texCoords = inTexCoords;
}
