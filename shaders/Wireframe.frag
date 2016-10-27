#version 410 core
in vec3 barycenter;

out vec4 color;

float edgeFactor()
{
    vec3 d = fwidth(barycenter);
    vec3 a3 = smoothstep(vec3(0.0), d, barycenter);
    return min(min(a3.x, a3.y), a3.z);
}

void main()
{
    color = vec4(1.0, 1.0, 1.0, (1.0-edgeFactor())*0.75);
}
