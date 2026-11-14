#define STACK_SIZE 160
#define GROUP_SIZE 64
#define NODE_LIST_SIZE 1024
#define LEAF_LIST_SIZE 512


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

    __local float4 nodeList[NODE_LIST_SIZE];    // accepted nodes: (cm, mass)
    __local int2   leafList[LEAF_LIST_SIZE];    // near leaves: (first body, count)
    __local int nodeCount;
    __local int leafCount;
    __local int walkDone;
    int stack[STACK_SIZE];
    int top = 0;
    stack[top++] = 0;


    int i = get_global_id(0);
    int lid = get_local_id(0);
    int g = get_group_id(0);
    float4 bi = (float4)(1.0f, 1.0f, 1.0f, 1.0f);
    bool valid = i < n;          
    if(valid) bi = posMass[i];
    float3 a = 0;

    float3 center;
    float3 boxHalf;

    if(lid == 0) {
        float3 lo = FLT_MAX;
        float3 hi = -FLT_MAX;
        for(int j = g*GROUP_SIZE; j < g*GROUP_SIZE+GROUP_SIZE && j < n ; ++j) {
            float3 pos = posMass[j].xyz;
            lo = fmin(lo, pos);
            hi = fmax(hi, pos);
        }
        center = (lo + hi) * 0.5f;
        boxHalf = 0.5f * (hi - lo);
    }

    while(true) {
        if(lid == 0) {
            nodeCount = 0;
            leafCount = 0;
            while(top > 0) {
                int idx = stack[--top];
                Node node = nodes[idx];
                if (node.meta.z == 0) continue;
                float3 boxDiff = fmax(fabs(center - node.cmMass.xyz) - boxHalf, 0.0f);
                float boxD2 = dot(boxDiff, boxDiff);
                float3 cellDiff = fmax(fabs(center - node.centerSize.xyz) - boxHalf, 0.0f);
                float cellD2 = dot(cellDiff, cellDiff);
                bool acceptBox = (node.centerSize.w * node.centerSize.w < boxD2 * theta2) && (cellD2 > node.centerSize.w * node.centerSize.w);
                if(node.meta.y == -1) {
                    if(leafCount + 1 > LEAF_LIST_SIZE) {
                        stack[top++] = idx;
                        break;
                    }
                    leafList[leafCount++] = (int2)(node.meta.x, node.meta.z);
                }
                else if(acceptBox) {
                    if(nodeCount + 1 > NODE_LIST_SIZE) {
                        stack[top++] = idx;
                        break;
                    }
                    nodeList[nodeCount++] = node.cmMass;
                }    
                else {
                    for(int c = 0; c < 8; ++c) {
                        stack[top++] = node.meta.y + c;
                    } 
                }
            }
            walkDone = top == 0;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        int nNodes = nodeCount;
        int nLeaves = leafCount;
        int done = walkDone;

        // accepted nodes: read from __local
        for(int k = 0; k < nNodes; ++k) {
            float4 bj = nodeList[k];
            float3 diff = bj.xyz - bi.xyz;
            float d2 = dot(diff, diff) + eps2;
            float inv_d = rsqrt(d2);
            float inv_D3 = inv_d * inv_d * inv_d;
            a += diff * (G * bj.w * inv_D3);
        }

        // near leaves: every lane reads the same posMass[j] at the same time (broadcast)
        for(int k = 0; k < nLeaves; ++k) {
            int2 leaf = leafList[k];
            for(int j = leaf.x; j < leaf.x + leaf.y; ++j) {
                float4 bj = posMass[j];
                float3 diff = bj.xyz - bi.xyz;
                float d2 = dot(diff, diff) + eps2;
                float inv_d = rsqrt(d2);
                float inv_D3 = inv_d * inv_d * inv_d;
                a += diff * (G * bj.w * inv_D3);
            }
        }
        
        barrier(CLK_LOCAL_MEM_FENCE);
        if(done) break;

    }

    if(valid) acc[i] = (float4)(a, 0.0f);
}
