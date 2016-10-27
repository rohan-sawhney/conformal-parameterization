#version 410 core
in vec3 position;
in vec3 normal;
in vec3 color;
in vec2 texCoords;

out vec4 fragColor;

layout (std140) uniform Light {
    vec3 position;
    vec3 color;
} light;

uniform vec3 viewPosition;

const vec3 BACKGROUND = vec3(0.1, 0.1, 0.1);
const vec3 ONE = vec3(1.0, 1.0, 1.0);

void main()
{
    vec3 lightDirection = normalize(light.position - position);
    vec3 viewDirection = normalize(viewPosition - position);
    vec3 halfwayDirection = normalize(lightDirection + viewDirection);
    
    float diffuse = max(dot(normal, lightDirection), 0.0);
    float specular = pow(max(dot(normal, halfwayDirection), 0.0), 12);
    float nv = max(dot(normal, viewDirection), 0.0);
    float fresnel = pow(sqrt(1.0 - nv*nv), 10);
    
    fragColor.rgb = (0.25 + 0.9*diffuse)*color + 0.2*specular*ONE + 0.75*fresnel*BACKGROUND;
    fragColor.a = 1.0;
}
