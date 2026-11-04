__kernel void directGravity(__global const float4* posMass,
                            __global float4*       acc,
                            const int   n,
                            const float G,
                            const float eps2)

{
    int i = get_global_id(0);
    if (i >= n) return;          

    
    float4 bi = posMass[i];

    float3 a = 0;

    for(int j = 0; j < n; ++j) {
        float4 bj = posMass[j];
        float3 d = bj.xyz - bi.xyz;
        float r2 = dot(d, d) + eps2;
        float inv_r = rsqrt(r2);
        float inv_r3 = inv_r*inv_r*inv_r;

        a += d * (G * bj.w * inv_r3);
    }

    acc[i] = (float4)(a, 0.0f);

}
