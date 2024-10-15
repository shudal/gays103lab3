using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class FVM : MonoBehaviour
{
	float dt 			= 0.003f;
    float mass 			= 1;
	float stiffness_0	= 20000.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;

	int[] 		Tet;
	int tet_number;			//The number of tetrahedra

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;
	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;

	//For Laplacian smoothing.
	Vector3[]   V_sum;
	int[]		V_num;

	SVD svd = new SVD();

    bool bDebugUseOnVolume = false;
    private float blendAlpha = 0.5f;

    private Vector3 floorPos = new Vector3(0, -3, 0);
    private Vector3 floorNormal = new Vector3(0, 1, 0);
    private float muN = 0.5f;
    private float muT = 0.5f;
    // Start is called before the first frame update
    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);
        	Tet = new int[tet_number*4];

    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}

        if (bDebugUseOnVolume)
        { 
            tet_number = 1;
            Tet = new int[tet_number * 4];
            Tet[0] = 0;
            Tet[1] = 1;
            Tet[2] = 2;
            Tet[3] = 3;

            number = 4;
            X = new Vector3[number];
            V = new Vector3[number];
            Force = new Vector3[number];
            X[0] = new Vector3(0, 0, 0);
            X[1] = new Vector3(1, 0, 0);
            X[2] = new Vector3(0, 1, 0);
            X[3] = new Vector3(0, 0, 1);
        }


        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

        //TODO: Need to allocate and assign inv_Dm
        inv_Dm = new Matrix4x4[tet_number];
        for (int i = 0; i < tet_number; i++)
        {
            inv_Dm[i] = Build_Edge_Matrix(i).inverse;
        }
    }
    Matrix4x4 Scale(Matrix4x4 m, float x)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int k = 0; k < 4; k++)
            {
                m[i, k] = m[i, k] * x;
            }
        }
        return m;
    }
    Matrix4x4 Add(Matrix4x4 m, Matrix4x4 n)
    {
        Matrix4x4 ans = m;
        for (int i = 0; i < 4; i++)
        {
            for (int k = 0; k < 4; k++)
            {
                ans[i, k] = m[i, k] + n[i, k];
            }
        }
        return ans;
    }
    Matrix4x4 BuildMatrixByThreeVector(Vector3 X10, Vector3 X20, Vector3 X30)
    {
        Matrix4x4 ret = Matrix4x4.zero;

        ret[0, 0] = X10.x;
        ret[1, 0] = X10.y;
        ret[2, 0] = X10.z;

        ret[0, 1] = X20.x;
        ret[1, 1] = X20.y;
        ret[2, 1] = X20.z;

        ret[0, 2] = X30.x;
        ret[1, 2] = X30.y;
        ret[2, 2] = X30.z;

        ret[3, 3] = 1f;
        return ret;
    }
    Matrix4x4 Build_Edge_Matrix(int tet)
    {
        Matrix4x4 ret = Matrix4x4.zero;
        //TODO: Need to build edge matrix here.
        if (tet >= 0 && tet < tet_number)
        {
            Vector3 X0 = X[Tet[tet * 4 + 0]];
            Vector3 X1 = X[Tet[tet * 4 + 1]];
            Vector3 X2 = X[Tet[tet * 4 + 2]];
            Vector3 X3 = X[Tet[tet * 4 + 3]];

            Vector3 X10 = X1 - X0;
            Vector3 X20 = X2 - X0;
            Vector3 X30 = X3 - X0;

            ret = BuildMatrixByThreeVector(X10, X20, X30);
        }
        return ret;
    }

    float trace(Matrix4x4 m)
    {
        float Res = 0;
        for (int i = 0; i < 4; i++)
        {
            Res += m[i, i];
        }
        return Res;
    }
    void _Update()
    {
        // Jump up.
        if (Input.GetKeyDown(KeyCode.Space))
        {
            for (int i = 0; i < number; i++)
                V[i].y += 0.2f;
        }

        for (int i = 0; i < number; i++)
        {
            //TODO: Add gravity to Force. 
            Force[i] = new Vector3(0, mass * -9.8f, 0);
        }

        for (int tet = 0; tet < tet_number; tet++)
        {
            //TODO: Deformation Gradient 
            int idx0 = Tet[tet * 4 + 0];
            int idx1 = Tet[tet * 4 + 1];
            int idx2 = Tet[tet * 4 + 2];
            int idx3 = Tet[tet * 4 + 3];


            if (bDebugUseOnVolume)
            {
                X[idx0] = new Vector3(-1f, 0f, 0f);
            }

            Vector3 X0 = X[idx0];
            Vector3 X1 = X[idx1];
            Vector3 X2 = X[idx2];
            Vector3 X3 = X[idx3];

            Vector3 X10 = X1 - X0;
            Vector3 X20 = X2 - X0;
            Vector3 X30 = X3 - X0;
            Matrix4x4 F = BuildMatrixByThreeVector(X10, X20, X30) * inv_Dm[tet];
            //TODO: Green Strain
            Matrix4x4 G = Scale(Add(F.transpose * F, Scale(Matrix4x4.identity, -1f)),0.5f);
            
            //TODO: Second PK Stress
            Matrix4x4 S = Add(Scale(G,2 * stiffness_1),Scale(Matrix4x4.identity, stiffness_0 * trace(G)));
            //TODO: Elastic Force
            Matrix4x4 f = Scale((F * S) * inv_Dm[tet].transpose, -(1f / 6f) * (1f / inv_Dm[tet].determinant));
            //f = Scale(f, 1f / f[3, 3]);

            Vector3 f1 = new Vector3(f[0, 0], f[1, 0], f[2, 0]);
            Vector3 f2 = new Vector3(f[0, 1], f[1, 1], f[2, 1]);
            Vector3 f3 = new Vector3(f[0, 2], f[1, 2], f[2, 2]);
            Vector3 f0 = -1f * (f1 + f2 + f3);

            Force[idx0] += f0;
            Force[idx1] += f1;
            Force[idx2] += f2;
            Force[idx3] += f3;
        }

        V_sum = new Vector3[V.Length];
        for (int i=0; i<V.Length; i++)
        {
            V_sum[i] = V[i];
            V_num[i] = 1;
        }
        for (int tet = 0; tet < tet_number; tet++)
        {
            //TODO: Deformation Gradient 
            int idx0 = Tet[tet * 4 + 0];
            int idx1 = Tet[tet * 4 + 1];
            int idx2 = Tet[tet * 4 + 2];
            int idx3 = Tet[tet * 4 + 3];
            Vector3 V0 = V[idx0];
            Vector3 V1 = V[idx1];
            Vector3 V2 = V[idx2];
            Vector3 V3 = V[idx3];

            V_sum[idx0] = V_sum[idx0] + V[idx1] + V[idx2] + V[idx3];
            V_num[idx0] += 3;

            V_sum[idx1] = V_sum[idx1] + V[idx0] + V[idx2] + V[idx3];
            V_num[idx1] += 3;

            V_sum[idx2] = V_sum[idx2] + V[idx1] + V[idx0] + V[idx3];
            V_num[idx2] += 3;

            V_sum[idx3] = V_sum[idx3] + V[idx1] + V[idx2] + V[idx0];
            V_num[idx3] += 3;
        }
        for (int i=0; i<V.Length; i++)
        {
            //V[i] = V_sum[i] / V_num[i];
            V[i] = blendAlpha * V[i] + (1 - blendAlpha) * V_sum[i] / V_num[i];
        }
        for (int i = 0; i < number; i++)
        {
            //TODO: Update X and V here.
            V[i] = V[i] * damp;
            V[i] += Force[i] * dt;
            X[i] = X[i] + V[i] * dt;
            //TODO: (Particle) collision with floor.
            /*
            Vector3 floorPos = GameObject.Find("Floor").transform.position;
            if (X[i].y < floorPos.y)
            {
                Vector3 posNew = X[i];
                posNew.y = floorPos.y;

                X[i] = posNew;
                V[i] = V[i] + (posNew - X[i]) * dt;
            }*/
            float signedDis = Vector3.Dot(X[i] - floorPos, floorNormal);
            if (signedDis < 0 && Vector3.Dot(V[i], floorNormal) < 0)
            {
                X[i] -= signedDis * floorNormal;
                Vector3 vN = Vector3.Dot(V[i], floorNormal) * floorNormal;
                Vector3 vT = V[i] - vN;
                float a = Math.Max(1 - muT * (1 + muN) * vN.magnitude / vT.magnitude, 0);
                V[i] = -muN * vN + a * vT;
            }
        }
    }

    // Update is called once per frame
    void Update()
    {
    	for(int l=0; l<10; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals ();
    }
}
