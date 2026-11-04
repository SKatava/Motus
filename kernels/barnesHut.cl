#define STACK_SIZE 160

typedef struct Node {
    float4   cmMass;         //4*4 bytes = 16
    float4   centerSize;     //4*4 bytes = 16 
    //meta.w = padding, .z = count, .y = firstChild, .x = first
    int4     meta;           //4*4 bytes = 16 => total = 48 bytes
} Node;


__kernel void barnesHutGravity(__global const float4* posMass,
                            __global float4*       acc,
                            __global const Node* nodes,
                            const int   n,
                            const float G,
                            const float eps2,
                            const float theta2)

{
    int i = get_global_id(0);
    if (i >= n) return;          
    float4 bi = posMass[i];
    float3 a = 0;

    int stack[STACK_SIZE];
    int top = 0;
    stack[top++] = 0;
     
    while(top > 0) {
        int idx = stack[--top];
        Node node = nodes[idx];
        if (node.meta.z == 0) continue;
        float3 diff = node.cmMass.xyz - bi.xyz;
        float d2 = dot(diff, diff) + eps2;
        float3 r = bi.xyz - node.centerSize.xyz; 
        if(node.meta.y == -1) {
            for(int j = node.meta.x; j < node.meta.x + node.meta.z; ++j) {
                if(j == i) continue;
                float4 bj = posMass[j];
                float3 diff = bj.xyz - bi.xyz;
                float d2 = dot(diff, diff) + eps2;
                float inv_d = rsqrt(d2);
                float inv_D3 = inv_d * inv_d * inv_d;
                a += diff * (G * bj.w * inv_D3);
            }
        }
        else if((node.centerSize.w * node.centerSize.w / d2 < theta2) && (dot(r,r) > node.centerSize.w * node.centerSize.w)) {
            float inv_d = rsqrt(d2);
            float inv_D3 = inv_d * inv_d * inv_d;

            a += diff * (G * node.cmMass.w * inv_D3);
        }    
        else {
            for(int c = 0; c < 8; ++c) {
                stack[top++] = node.meta.y + c;
            } 
        }
    }

    acc[i] = (float4)(a, 0.0f);
}
