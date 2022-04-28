//Test Assimp Application
//Import Model Data
#include <algorithm>
#include <iostream>
#include <filesystem>


#include "DataFileContainer.h"
#include "TransformNode.h"
#include "AnimatedMesh.h"

void RemoveCharsFromString(std::string& str, char charToRemove)
{
	str.erase(std::ranges::remove(str, charToRemove).begin(), str.end());
}

void printEndLine(float lines)
{
	for (int i = 0; i < lines; i++)
	{
		std::cout << std::endl;
	}
}

void printValueToConsole(std::string valHeader, int val)
{
	std::cout << "  " << valHeader << ": " << val << std::endl;
}

void printValueToConsole(std::string valHeader, bool val)
{
	std::cout << "  " << valHeader << ": " << val << std::endl;
}

void printValueToConsole(std::string valHeader, float val)
{
	std::cout << "  " << valHeader << ": " << val << std::endl;
}

void printValueToConsole(std::string valHeader, std::string val)
{
	std::cout << "  " << valHeader << ": " << val << std::endl;
}

void printValueToConsole(std::string valHeader, Vertex val)
{
	std::cout << "  " << valHeader << " Position X: " << val.Position.x << std::endl;
	std::cout << "  " << valHeader << " Position Y: " << val.Position.y << std::endl;
	std::cout << "  " << valHeader << " Position Z: " << val.Position.z << std::endl;

	std::cout << "  " << valHeader << " Normal X: " << val.Normal.x << std::endl;
	std::cout << "  " << valHeader << " Normal Y: " << val.Normal.y << std::endl;
	std::cout << "  " << valHeader << " Normal Z: " << val.Normal.z << std::endl;

	std::cout << "  " << valHeader << " TexCoord X: " << val.TexCoords.x << std::endl;
	std::cout << "  " << valHeader << " TexCoord Y: " << val.TexCoords.y << std::endl;
}
void printVerticesFromVectorToConsole(std::string valHeader, std::vector<float>& val)
{
	int vertexcount = val.size() / 8;
	Vertex vertices;

	for (int k = 0; k < val.size(); k++)
	{
		std::cout << " " << val[k] << " " << val[k + 1] << " " << val[k + 2] << " " << val[k + 3]
			<< " " << val[k + 4] << " " << val[k + 5] << " " << val[k + 6] << " " << val[k + 7] << std::endl;
		std::cout << std::endl;

		k += 7;
	}
}


struct KeySRT
{
	KeySRT(const KeySRT& other)
		: Timestamp(other.Timestamp),
		  mGlobalTransform(other.mGlobalTransform)
	{
	}

	KeySRT(KeySRT&& other) noexcept
		: Timestamp(other.Timestamp),
		  mGlobalTransform(std::move(other.mGlobalTransform))
	{
	}

	KeySRT& operator=(const KeySRT& other)
	{
		if (this == &other)
			return *this;
		Timestamp = other.Timestamp;
		mGlobalTransform = other.mGlobalTransform;
		return *this;
	}

	KeySRT& operator=(KeySRT&& other) noexcept
	{
		if (this == &other)
			return *this;
		Timestamp = other.Timestamp;
		mGlobalTransform = std::move(other.mGlobalTransform);
		return *this;
	}

	KeySRT() = default;

	float Timestamp;
	float PositionX;
	float PositionY;
	float PositionZ;
	float RotationX;
	float RotationY;
	float RotationZ;
	float RotationW;
	float ScalingX;
	float ScalingY;
	float ScalingZ;

	FbxAMatrix mGlobalTransform;
};

struct FBXBoneAnimationKeyframes
{
	FBXBoneAnimationKeyframes(const FBXBoneAnimationKeyframes& other)
		: boneKeyframes(other.boneKeyframes),
		  BoneName(other.BoneName),
		  BoneId(other.BoneId)
	{
	}
	FBXBoneAnimationKeyframes(FBXBoneAnimationKeyframes&& other) noexcept
		: boneKeyframes(std::move(other.boneKeyframes)),
		  BoneName(std::move(other.BoneName)),
		  BoneId(other.BoneId)
	{
	}
	FBXBoneAnimationKeyframes& operator=(const FBXBoneAnimationKeyframes& other)
	{
		if (this == &other)
			return *this;
		boneKeyframes = other.boneKeyframes;
		BoneName = other.BoneName;
		BoneId = other.BoneId;
		return *this;
	}
	FBXBoneAnimationKeyframes& operator=(FBXBoneAnimationKeyframes&& other) noexcept
	{
		if (this == &other)
			return *this;
		boneKeyframes = std::move(other.boneKeyframes);
		BoneName = std::move(other.BoneName);
		BoneId = other.BoneId;
		return *this;
	}
	FBXBoneAnimationKeyframes() = default;


	std::vector<KeySRT> boneKeyframes;
	std::string BoneName;
	int BoneId;
};

struct Joint
{
	std::string mName;
	int mParentIndex;
	int BoneId;
	FbxAMatrix mGlobalBindposeInverse;
	FbxNode* mNode;

	Joint() :
		mNode(nullptr)
	{
		mGlobalBindposeInverse.SetIdentity();
		mParentIndex = -1;
	}

	~Joint()
	{
	}
};

struct FBXAnimation
{
	FBXAnimation(const FBXAnimation& other)
		: AnimationKeyframes(other.AnimationKeyframes),
		  AnimationName(other.AnimationName),
		  AnimationLength(other.AnimationLength)
	{
	}

	FBXAnimation(FBXAnimation&& other) noexcept
		: AnimationKeyframes(other.AnimationKeyframes),
		  AnimationName(std::move(other.AnimationName)),
		  AnimationLength(other.AnimationLength)
	{
	}

	FBXAnimation& operator=(const FBXAnimation& other)
	{
		if (this == &other)
			return *this;
		AnimationKeyframes = other.AnimationKeyframes;
		AnimationName = other.AnimationName;
		AnimationLength = other.AnimationLength;
		return *this;
	}

	FBXAnimation& operator=(FBXAnimation&& other) noexcept
	{
		if (this == &other)
			return *this;
		AnimationKeyframes = other.AnimationKeyframes;
		AnimationName = std::move(other.AnimationName);
		AnimationLength = other.AnimationLength;
		return *this;
	}

	std::vector<FBXBoneAnimationKeyframes> AnimationKeyframes;
	std::string AnimationName;
	float AnimationLength;
	FBXAnimation()
	{
	}

	~FBXAnimation()
	{
	}
};

struct Skeleton
{
	std::vector<Joint> mJoints;
};


struct FBXSceneNode
{
	FBXSceneNode(const FBXSceneNode& other)
		: ParentNodeID(other.ParentNodeID),
		  NodeID(other.NodeID),
		  NodeTransform(other.NodeTransform),
		  NodeName(other.NodeName),
		  NumberOfChilds(other.NumberOfChilds),
		  NumberOfMeshes(other.NumberOfMeshes),
		  AnimatedNodeMeshes(other.AnimatedNodeMeshes),
		  ChildNodes(other.ChildNodes),
		  mMaterialLookUp(other.mMaterialLookUp)
	{
	}

	FBXSceneNode(FBXSceneNode&& other) noexcept
		: ParentNodeID(other.ParentNodeID),
		  NodeID(other.NodeID),
		  NodeTransform(std::move(other.NodeTransform)),
		  NodeName(std::move(other.NodeName)),
		  NumberOfChilds(other.NumberOfChilds),
		  NumberOfMeshes(other.NumberOfMeshes),
		  AnimatedNodeMeshes(std::move(other.AnimatedNodeMeshes)),
		  ChildNodes(std::move(other.ChildNodes)),
		  mMaterialLookUp(std::move(other.mMaterialLookUp))
	{
	}

	FBXSceneNode& operator=(const FBXSceneNode& other)
	{
		if (this == &other)
			return *this;
		ParentNodeID = other.ParentNodeID;
		NodeID = other.NodeID;
		NodeTransform = other.NodeTransform;
		NodeName = other.NodeName;
		NumberOfChilds = other.NumberOfChilds;
		NumberOfMeshes = other.NumberOfMeshes;
		AnimatedNodeMeshes = other.AnimatedNodeMeshes;
		ChildNodes = other.ChildNodes;
		mMaterialLookUp = other.mMaterialLookUp;
		return *this;
	}

	FBXSceneNode& operator=(FBXSceneNode&& other) noexcept
	{
		if (this == &other)
			return *this;
		ParentNodeID = other.ParentNodeID;
		NodeID = other.NodeID;
		NodeTransform = std::move(other.NodeTransform);
		NodeName = std::move(other.NodeName);
		NumberOfChilds = other.NumberOfChilds;
		NumberOfMeshes = other.NumberOfMeshes;
		AnimatedNodeMeshes = std::move(other.AnimatedNodeMeshes);
		ChildNodes = std::move(other.ChildNodes);
		mMaterialLookUp = std::move(other.mMaterialLookUp);
		return *this;
	}

private:
	static void ComputeRandomColor(HACD::Material& mat)
	{
		mat.m_diffuseColor.X() = mat.m_diffuseColor.Y() = mat.m_diffuseColor.Z() = 0.0f;
		while (mat.m_diffuseColor.X() == mat.m_diffuseColor.Y() || mat.m_diffuseColor.Z() == mat.m_diffuseColor.Y() ||
			mat.m_diffuseColor.Z() == mat.m_diffuseColor.X())
		{
			mat.m_diffuseColor.X() = (rand() % 100) / 100.0f;
			mat.m_diffuseColor.Y() = (rand() % 100) / 100.0f;
			mat.m_diffuseColor.Z() = (rand() % 100) / 100.0f;
		}
	}

	static bool SaveVRML2(std::ofstream& fout, const double* const& points, const int* const& triangles,
	                      const unsigned& nPoints, const unsigned& nTriangles, const HACD::Material& material)
	{
		if (fout.is_open())
		{
			fout.setf(std::ios::fixed, std::ios::floatfield);
			fout.setf(std::ios::showpoint);
			fout.precision(6);
			size_t nV = nPoints * 3;
			size_t nT = nTriangles * 3;
			fout << "#VRML V2.0 utf8" << std::endl;
			fout << "" << std::endl;
			fout << "# Vertices: " << nPoints << std::endl;
			fout << "# Triangles: " << nTriangles << std::endl;
			fout << "" << std::endl;
			fout << "Group {" << std::endl;
			fout << "    children [" << std::endl;
			fout << "        Shape {" << std::endl;
			fout << "            appearance Appearance {" << std::endl;
			fout << "                material Material {" << std::endl;
			fout << "                    diffuseColor " << material.m_diffuseColor.X() << " "
				<< material.m_diffuseColor.Y() << " "
				<< material.m_diffuseColor.Z() << std::endl;
			fout << "                    ambientIntensity " << material.m_ambientIntensity << std::endl;
			fout << "                    specularColor " << material.m_specularColor.X() << " "
				<< material.m_specularColor.Y() << " "
				<< material.m_specularColor.Z() << std::endl;
			fout << "                    emissiveColor " << material.m_emissiveColor.X() << " "
				<< material.m_emissiveColor.Y() << " "
				<< material.m_emissiveColor.Z() << std::endl;
			fout << "                    shininess " << material.m_shininess << std::endl;
			fout << "                    transparency " << material.m_transparency << std::endl;
			fout << "                }" << std::endl;
			fout << "            }" << std::endl;
			fout << "            geometry IndexedFaceSet {" << std::endl;
			fout << "                ccw TRUE" << std::endl;
			fout << "                solid TRUE" << std::endl;
			fout << "                convex TRUE" << std::endl;
			if (nV > 0)
			{
				fout << "                coord DEF co Coordinate {" << std::endl;
				fout << "                    point [" << std::endl;
				for (size_t v = 0; v < nV; v += 3)
				{
					fout << "                        " << points[v + 0] << " "
						<< points[v + 1] << " "
						<< points[v + 2] << "," << std::endl;
				}
				fout << "                    ]" << std::endl;
				fout << "                }" << std::endl;
			}
			if (nT > 0)
			{
				fout << "                coordIndex [ " << std::endl;
				for (size_t f = 0; f < nT; f += 3)
				{
					fout << "                        " << triangles[f + 0] << ", "
						<< triangles[f + 1] << ", "
						<< triangles[f + 2] << ", -1," << std::endl;
				}
				fout << "                ]" << std::endl;
			}
			fout << "            }" << std::endl;
			fout << "        }" << std::endl;
			fout << "    ]" << std::endl;
			fout << "}" << std::endl;
			return true;
		}
		else
		{
			return false;
		}
	}

	struct DecompResults
	{
		btCompoundShape* CompoundShape = nullptr;
		btAlignedObjectArray<btConvexShape*> m_convexShapes = {};
		btAlignedObjectArray<btTriangleMesh*> m_trimeshes = {};
	};

	static DecompResults* Decomp(std::vector<int> indicies, std::vector<float> verticePositions)
	{
		return nullptr;
		/*//
		//	Setup Indices
		const uint32_t nTriangles = indicies.size();
		std::vector<uint32_t> Triangles;
		for (uint32_t i = 0; i < nTriangles; i++)
		{
			Triangles.push_back(indicies[i]);
		}
		//
		//	Setup Points (3 Points is 1 Vertex)
		const uint32_t nPoints = verticePositions.size();
		std::vector<float> Points;
		for (uint32_t i = 0; i < nPoints; i++)
		{
			Points.push_back(verticePositions[i]);
			i++;
			Points.push_back(verticePositions[i]);
			i++;
			Points.push_back(verticePositions[i]);
		}
		//
		//	Setup VHACD Parameters and create its interface
		VHACD::IVHACD::Parameters params;
		VHACD::IVHACD* interfaceVHACD = VHACD::CreateVHACD();
		VHACD::IVHACD::ConvexHull Hull;
		//
		//	Compute approximate convex decomposition
		//printf("Compute V-HACD: Points %i Triangles %i\n", Points.size(), Triangles.size());
		bool res = interfaceVHACD->Compute(Points.data(), (uint32_t)(Points.size() / 3),
		                                   Triangles.data(), (uint32_t)(Triangles.size() / 3), params);
		//
		//	Get the number of convex hulls
		unsigned int nConvexHulls = interfaceVHACD->GetNConvexHulls();
		//printf("V-HACD Done: Hull Count %i\n", nConvexHulls);
		//
		//	Create a new DecompResults structure
		DecompResults* Results = new DecompResults;
		//
		//	Create a new Compound Shape for this decomposition
		Results->CompoundShape = new btCompoundShape();
		//
		//	Iterate through each convex hull and fill results,
		for (unsigned int h = 0; h < nConvexHulls; ++h)
		{
			//printf("\tHull: %i\n", h);
			//printf("\t\tPoints: %i\n", Hull.m_points);
			//printf("\t\tTriangles: %i\n", Hull.m_triangles);
			//printf("\t\tVertices: %i\n", vertices.size());
			//
			//	Fill 'Hull' for each individual convex hull
			interfaceVHACD->GetConvexHull(h, Hull);
			//
			//	Create a new Triangle Mesh for this hull
			btTriangleMesh* trimesh = new btTriangleMesh();
			Results->m_trimeshes.push_back(trimesh);
			//
			//	Grab the hulls center position
			const btVector3 centroid(Hull.m_center[0], Hull.m_center[1], Hull.m_center[2]);
			printf("Hull Center %f %f %f\n", Hull.m_center[0], Hull.m_center[1], Hull.m_center[2]);
			//
			//	Iterate through this hulls triangles
			for (unsigned int i = 0; i < Hull.m_nTriangles; i++)
			{
				//
				//	Calculate indices
				const unsigned int index0 = Hull.m_triangles[i * 3];
				const unsigned int index1 = Hull.m_triangles[i * 3 + 1];
				const unsigned int index2 = Hull.m_triangles[i * 3 + 2];
				//
				//	Calculate vertices
				const btVector3 vertex0(Hull.m_points[index0 * 3], Hull.m_points[index0 * 3 + 1],
				                        Hull.m_points[index0 * 3 + 2]);
				const btVector3 vertex1(Hull.m_points[index1 * 3], Hull.m_points[index1 * 3 + 1],
				                        Hull.m_points[index1 * 3 + 2]);
				const btVector3 vertex2(Hull.m_points[index2 * 3], Hull.m_points[index2 * 3 + 1],
				                        Hull.m_points[index2 * 3 + 2]);
				//
				//	Add this triangle into our Triangle Mesh
				trimesh->addTriangle(vertex0 - centroid, vertex1 - centroid, vertex2 - centroid);
			}
			//
			//	Create a new ConvexShape from this hulls Triangle Mesh
			btConvexShape* convexShape = new btConvexTriangleMeshShape(trimesh);
			Results->m_convexShapes.push_back(convexShape);
			//
			//	Create a transform using centroid as origin
			btTransform trans;
			trans.setIdentity();
			trans.setOrigin(centroid);
			//
			//	Add this ConvexShape to our CompoundShape
			Results->CompoundShape->addChildShape(trans, convexShape);
		}
		std::ofstream foutCH("Test.wrl");

		VHACD::IVHACD::ConvexHull ch;
		if (foutCH.is_open())
		{
			HACD::Material mat;
			for (unsigned int p = 0; p < nConvexHulls; ++p)
			{
				interfaceVHACD->GetConvexHull(p, ch);
				ComputeRandomColor(mat);
				SaveVRML2(foutCH, ch.m_points, (const int*)ch.m_triangles, ch.m_nPoints, ch.m_nTriangles, mat);
			}
			foutCH.close();
		}

		//
		// release memory
		interfaceVHACD->Clean();
		interfaceVHACD->Release();
		//
		//	Return our DecompResults containing the CompoundShape full of Convexically Decomposed Convex Shapes
		return Results;*/
	}

	void GetSubsetData(FbxMesh* pMesh, std::vector<TriangleVertexIndicies>& OutIndicies, std::vector<SubsetSurfaceData>& OutSubsets)
	{
		FbxLayerElementArrayTemplate<int>* materialIndices;
		FbxGeometryElement::EMappingMode materialMappingMode = FbxGeometryElement::eNone;
		FbxMesh* currMesh = pMesh;

		if (currMesh->GetElementMaterial())
		{
			materialIndices = &(currMesh->GetElementMaterial()->GetIndexArray());
			materialMappingMode = currMesh->GetElementMaterial()->GetMappingMode();

			if (materialIndices)
			{
				switch (materialMappingMode)
				{
				case FbxGeometryElement::eByPolygon:
				{
					if (materialIndices->GetCount() == OutIndicies.size())
					{
						int currentMaterialIndex = -1;
						int currentMaterialStart = 0;
						for (unsigned int i = 0; i < OutIndicies.size(); ++i)
						{
							unsigned int materialIndex = materialIndices->GetAt(i);
							if (materialIndex != currentMaterialIndex)
							{
								OutSubsets.push_back(SubsetSurfaceData());
								OutSubsets.back().SurfaceMaterial.materialIndex = currentMaterialIndex;
								OutSubsets.back().Start = currentMaterialStart * 3;
								OutSubsets.back().DrawAmount = (i - currentMaterialStart) * 3;
								currentMaterialIndex = materialIndex;
								currentMaterialStart = i;
							}
						}
					}
				}
				break;

				case FbxGeometryElement::eAllSame:
				{
					unsigned int materialIndex = materialIndices->GetAt(0);
					OutSubsets.push_back(SubsetSurfaceData());
					OutSubsets.back().SurfaceMaterial.materialIndex = materialIndex;
					OutSubsets.back().Start = 0;
					OutSubsets.back().DrawAmount = OutIndicies.size() * 3;
				}
				break;

				default:
					throw std::exception("Invalid mapping mode for material\n");
				}
			}
		}
	}

	static unsigned int FindJointIndexUsingName(const std::string& inJointName, Skeleton& skeleton)
	{
		for (unsigned int i = 0; i < skeleton.mJoints.size(); ++i)
		{
			if (skeleton.mJoints[i].mName == inJointName)
			{
				return skeleton.mJoints[i].BoneId;
			}
		}

		throw std::exception("Skeleton information in FBX file is corrupted.");
	}

	static void GetAnimationsBoneIdsAndBoneWeights(FbxNode* inNode, std::vector<AnimatedVertex>& OutVertices, Skeleton& skeleton, std::vector<FBXAnimation>& animations)
	{
		FbxMesh* currMesh = inNode->GetMesh();
		unsigned int numOfDeformers = currMesh->GetDeformerCount();
		// This geometry transform is something I cannot understand
		// I think it is from MotionBuilder
		// If you are using Maya for your models, 99% this is just an
		// identity matrix
		// But I am taking it into account anyways......
		//FbxAMatrix geometryTransform = Utilities::GetGeometryTransformation(inNode);

		const FbxVector4 lT = inNode->GetGeometricTranslation(FbxNode::eSourcePivot);
		const FbxVector4 lR = inNode->GetGeometricRotation(FbxNode::eSourcePivot);
		const FbxVector4 lS = inNode->GetGeometricScaling(FbxNode::eSourcePivot);
		FbxAMatrix geometryTransform = FbxAMatrix(lT, lR, lS);
		// A deformer is a FBX thing, which contains some clusters
		// A cluster contains a link, which is basically a joint
		// Normally, there is only one deformer in a mesh
		for (unsigned int deformerIndex = 0; deformerIndex < numOfDeformers; ++deformerIndex)
		{
			// There are many types of deformers in Maya,
			// We are using only skins, so we see if this is a skin
			FbxSkin* currSkin = reinterpret_cast<FbxSkin*>(currMesh->GetDeformer(deformerIndex, FbxDeformer::eSkin));
			if (!currSkin)
			{
				continue;
			}

			unsigned int numOfClusters = currSkin->GetClusterCount();
			for (unsigned int clusterIndex = 0; clusterIndex < numOfClusters; ++clusterIndex)
			{
				FbxCluster* currCluster = currSkin->GetCluster(clusterIndex);
				std::string currJointName = currCluster->GetLink()->GetName();
				unsigned int currJointIndex = FindJointIndexUsingName(currJointName, skeleton);
				FbxAMatrix transformMatrix;
				FbxAMatrix transformLinkMatrix;
				FbxAMatrix globalBindposeInverseMatrix;

				currCluster->GetTransformMatrix(transformMatrix); // The transformation of the mesh at binding time
				currCluster->GetTransformLinkMatrix(transformLinkMatrix);


				globalBindposeInverseMatrix = transformLinkMatrix * transformMatrix * geometryTransform;
				// Update the information in mSkeleton 
				skeleton.mJoints[currJointIndex].mGlobalBindposeInverse = globalBindposeInverseMatrix.Inverse();
				skeleton.mJoints[currJointIndex].mNode = currCluster->GetLink();

				// Associate each joint with the control points it affects
				unsigned int numOfIndices = currCluster->GetControlPointIndicesCount();
				for (unsigned int i = 0; i < numOfIndices; ++i)
				{
					if (OutVertices[currCluster->GetControlPointIndices()[i]].BoneIDBoneWeight.size() < 4)
					{
						BoneIdBoneWeight temp;
						temp.BoneId = currJointIndex;
						temp.BoneWeight = currCluster->GetControlPointWeights()[i];
						OutVertices[currCluster->GetControlPointIndices()[i]].BoneIDBoneWeight.push_back(temp);
					}

				}

				// Get animation information
				// Now only supports one take
				for(int animationStackIndex = 0; animationStackIndex < inNode->GetScene()->GetSrcObjectCount<FbxAnimStack>(); animationStackIndex++)
				{
					FBXAnimation& currentAnimation = animations[animationStackIndex];
					FbxAnimStack* currAnimStack = inNode->GetScene()->GetSrcObject<FbxAnimStack>(animationStackIndex);
					inNode->GetScene()->SetCurrentAnimationStack(currAnimStack);
					FbxAnimEvaluator* animEvaluator = currCluster->GetLink()->GetAnimationEvaluator();
					FbxString animStackName = currAnimStack->GetName();
					FbxTakeInfo* takeInfo = inNode->GetScene()->GetTakeInfo(animStackName);
					FbxTime start = takeInfo->mLocalTimeSpan.GetStart();
					FbxTime end = takeInfo->mLocalTimeSpan.GetStop();
					currentAnimation.AnimationKeyframes.push_back(FBXBoneAnimationKeyframes());
					FBXBoneAnimationKeyframes& currAnim = currentAnimation.AnimationKeyframes.back();
					currAnim.BoneId = skeleton.mJoints[currJointIndex].BoneId;
					currAnim.BoneName = skeleton.mJoints[currJointIndex].mName;
					FbxLongLong FrameCount = end.GetFrameCount(FbxTime::eFrames24);
					for (FbxLongLong i = start.GetFrameCount(FbxTime::eFrames24); i <= FrameCount
						; ++i)
					{
						FbxTime currTime;
						currTime.SetFrame(i, FbxTime::eFrames24);
						currAnim.boneKeyframes.push_back(KeySRT());
						currAnim.boneKeyframes.back().Timestamp = currTime.GetMilliSeconds();
						FbxAMatrix currentTransformOffset = inNode->EvaluateGlobalTransform(currTime) * geometryTransform;
						currAnim.boneKeyframes.back().mGlobalTransform =currentTransformOffset.Inverse() * currCluster->GetLink()->EvaluateGlobalTransform(currTime);
					}
					inNode->GetScene()->SetCurrentAnimationStack(nullptr);
				}

			
			}
		}

		for (auto itr = OutVertices.begin(); itr != OutVertices.end(); ++itr)
		{
			while (itr->BoneIDBoneWeight.size() < 4)
			{
				BoneIdBoneWeight temp;
				temp.BoneId = -1;
				temp.BoneWeight = 0.0f;
				itr->BoneIDBoneWeight.push_back(temp);
			}
		}
	}
	void ProcessMaterials(FbxNode* inNode)
	{
		unsigned int materialCount = inNode->GetMaterialCount();

		for (unsigned int i = 0; i < materialCount; ++i)
		{
			FbxSurfaceMaterial* surfaceMaterial = inNode->GetMaterial(i);
			ProcessMaterialTexture(surfaceMaterial, mMaterialLookUp[i]);
		}
	}
	void ProcessMaterialTexture(FbxSurfaceMaterial* inMaterial, SurfaceMaterialData& ioMaterial)
	{
		unsigned int textureIndex = 0;
		FbxProperty property;
		ioMaterial.matName = inMaterial->GetName();
		FBXSDK_FOR_EACH_TEXTURE(textureIndex)
		{
			property = inMaterial->FindProperty(FbxLayerElement::sTextureChannelNames[textureIndex]);
			if (property.IsValid())
			{
				unsigned int textureCount = property.GetSrcObjectCount<FbxTexture>();
				for (unsigned int i = 0; i < textureCount; ++i)
				{
					FbxLayeredTexture* layeredTexture = property.GetSrcObject<FbxLayeredTexture>(i);
					if (layeredTexture)
					{
						throw std::exception("Layered Texture is currently unsupported\n");
					}
					else
					{
						FbxTexture* texture = property.GetSrcObject<FbxTexture>(i);
						if (texture)
						{
							std::string textureType = property.GetNameAsCStr();
							FbxFileTexture* fileTexture = FbxCast<FbxFileTexture>(texture);

							if (fileTexture)
							{
								if (textureType == "DiffuseColor")
								{
									ioMaterial.AlbedoMapFilename = fileTexture->GetFileName();
								}
							}
						}
					}
				}
			}
		}
	}

public:
	FBXSceneNode(std::string nodeName);
	FBXSceneNode(std::string nodeName, uint64_t parentNodeId, TransformNode* parentTransform, uint64_t nodeId,
	             DirectX::FXMVECTOR position, DirectX::FXMVECTOR quaternionRotation, DirectX::FXMVECTOR scale);

	FBXSceneNode(std::string nodeName, uint64_t parentNodeId, TransformNode* parentTransform, uint64_t nodeId);

	uint64_t ParentNodeID;
	uint64_t NodeID;
	TransformNode NodeTransform;
	std::string NodeName;
	int NumberOfChilds{};
	int NumberOfMeshes{};
	std::vector<AnimatedMesh> AnimatedNodeMeshes{};
	std::vector<FBXSceneNode> ChildNodes{};
	std::unordered_map<unsigned int, SurfaceMaterialData> mMaterialLookUp;

	static void ReadNodeSkeletonData(FbxNode* pNode, int pDepth, int myIndex, int inParentIndex, Skeleton& skeleton)
	{
		if (pNode->GetNodeAttribute() && pNode->GetNodeAttribute()->GetAttributeType() && pNode->GetNodeAttribute()->
			GetAttributeType() == FbxNodeAttribute::eSkeleton)
		{
			Joint currJoint;
			currJoint.mName = pNode->GetName();
			RemoveCharsFromString(currJoint.mName, '.');
			currJoint.BoneId = myIndex;
			currJoint.mParentIndex = inParentIndex;
			skeleton.mJoints.push_back(currJoint);
		}

		for (int i = 0; i < pNode->GetChildCount(); i++)
		{
			ReadNodeSkeletonData(pNode->GetChild(i), pDepth + 1, skeleton.mJoints.size(), myIndex, skeleton);
		}
	}

	

	void ReadNodeMeshData(FbxNode* pNode, TransformNode* parent, int pDepth, Skeleton& skeleton, bool convertToLeftHanded, std::vector<FBXAnimation>& animations)
	{
		if (pNode->GetNodeAttribute())
		{
			const FbxVector4 lT = pNode->GetGeometricTranslation(FbxNode::eSourcePivot);
			const FbxVector4 lR = pNode->GetGeometricRotation(FbxNode::eSourcePivot);
			const FbxVector4 lS = pNode->GetGeometricScaling(FbxNode::eSourcePivot);

			DirectX::XMFLOAT3 pos;
			DirectX::XMFLOAT4 rotation;
			DirectX::XMFLOAT3 scale;

			pos.x = static_cast<float>(lT.mData[0]);
			pos.y = static_cast<float>(lT.mData[1]);
			pos.z = static_cast<float>(lT.mData[2]);

			rotation.x = static_cast<float>(lR.mData[0]);
			rotation.y = static_cast<float>(lR.mData[1]);
			rotation.z = static_cast<float>(lR.mData[2]);
			rotation.w = static_cast<float>(lR.mData[3]);

			scale.x = static_cast<float>(lS.mData[0]);
			scale.y = static_cast<float>(lS.mData[1]);
			scale.z = static_cast<float>(lS.mData[2]);
			std::string nodeName = pNode->GetName();
			RemoveCharsFromString(nodeName, '.');
			if (pNode->GetParent())
			{
				this->ParentNodeID = pNode->GetParent()->GetUniqueID();
				this->NodeTransform = TransformNode(nodeName, parent, DirectX::XMLoadFloat3(&pos),
				                                    DirectX::XMLoadFloat4(&rotation), DirectX::XMLoadFloat3(&scale));
			}
			else
			{
				this->ParentNodeID = -1;
				this->NodeTransform = TransformNode(nodeName, nullptr, DirectX::XMLoadFloat3(&pos),
				                                    DirectX::XMLoadFloat4(&rotation), DirectX::XMLoadFloat3(&scale));
			}

			this->NodeID = pNode->GetUniqueID();
			this->NodeName = nodeName;
			this->NumberOfChilds = pNode->GetChildCount();
			this->NumberOfMeshes = 0;

			FbxString lString;
			int i;

			for (i = 0; i < pDepth; i++)
			{
				lString += "     ";
			}

			lString += pNode->GetName();
			lString += "\n";

			FBXSDK_printf(lString.Buffer());

			for (int nodeAttributeIndex = 0; nodeAttributeIndex < 1; nodeAttributeIndex++)
			{
				FbxNodeAttribute* node_attribute = pNode->GetNodeAttributeByIndex(nodeAttributeIndex);

				if (node_attribute->GetAttributeType() == FbxNodeAttribute::EType::eMesh)
				{
					this->NumberOfMeshes = 1;
					FbxGeometryConverter lConverter(pNode->GetFbxManager());
					node_attribute = lConverter.Triangulate(node_attribute, true);

					auto* Mesh = static_cast<FbxMesh*>(node_attribute);

					std::vector<SubsetSurfaceData> SubsetData = {};
					std::vector<AnimatedVertex> OutVertices = {};
					std::vector<AnimatedVertex> TempVertices = {};
					std::vector<TriangleVertexIndicies> OutIndicies = {};
					int counter = 0;
					FbxVector4* Vertices = Mesh->GetControlPoints();
					for (int j = 0; j < Mesh->GetControlPointsCount(); j++)
					{
						AnimatedVertex NewVertex{};
						NewVertex.Position.x = (float)Vertices[j].mData[0];
						NewVertex.Position.y = (float)Vertices[j].mData[1];
						NewVertex.Position.z = (float)Vertices[j].mData[2];
						TempVertices.push_back(NewVertex);
					}

					GetAnimationsBoneIdsAndBoneWeights(pNode, TempVertices, skeleton, animations);
					for (int p = 0; p < Mesh->GetPolygonCount(); p++)
					{
						if(convertToLeftHanded)
						{
							int indexCounter = 0;
							OutIndicies.push_back(TriangleVertexIndicies());
							//Loop through the three vertices within each polygon
							for (int v = Mesh->GetPolygonSize(p) - 1; v >= 0; v--)
							{

								const int controlPointIndex = Mesh->GetPolygonVertex(p, v);

								FbxVector4 controlPointPosition = Mesh->GetControlPointAt(controlPointIndex);

								AnimatedVertex NewVertex{};
								NewVertex.Position.x = (float)controlPointPosition.mData[0];
								NewVertex.Position.y = (float)controlPointPosition.mData[1];
								NewVertex.Position.z = (float)controlPointPosition.mData[2] * -1.0f;

								FbxStringList lUVSetNameList;
								Mesh->GetUVSetNames(lUVSetNameList);

								FbxVector2 tex;
								bool isMapped;
								Mesh->GetPolygonVertexUV(p, v, lUVSetNameList.GetStringAt(0), tex, isMapped);
								NewVertex.TexCoords.x = (float)tex.mData[0];
								NewVertex.TexCoords.y = 1.0f - (float)tex.mData[1];

								FbxVector4 norm;
								Mesh->GetPolygonVertexNormal(p, v, norm);
								NewVertex.Normal.x = (float)norm.mData[0];
								NewVertex.Normal.y = (float)norm.mData[1];
								NewVertex.Normal.z = (float)norm.mData[2] * -1.0f;

								NewVertex.BoneIDBoneWeight = TempVertices[controlPointIndex].BoneIDBoneWeight;

								OutVertices.push_back(NewVertex);
								OutIndicies.back().Indices[indexCounter] = counter;
								indexCounter++;
								counter++;
							}
						}
						else
						{
							OutIndicies.push_back(TriangleVertexIndicies());
							//Loop through the three vertices within each polygon
							for (int v = 0; v < Mesh->GetPolygonSize(p); v++)
							{

								const int controlPointIndex = Mesh->GetPolygonVertex(p, v);

								FbxVector4 controlPointPosition = Mesh->GetControlPointAt(controlPointIndex);

								AnimatedVertex NewVertex{};
								NewVertex.Position.x = (float)controlPointPosition.mData[0];
								NewVertex.Position.y = (float)controlPointPosition.mData[1];
								NewVertex.Position.z = (float)controlPointPosition.mData[2];

								FbxStringList lUVSetNameList;
								Mesh->GetUVSetNames(lUVSetNameList);

								FbxVector2 tex;
								bool isMapped;
								Mesh->GetPolygonVertexUV(p, v, lUVSetNameList.GetStringAt(0), tex, isMapped);
								NewVertex.TexCoords.x = (float)tex.mData[0];
								NewVertex.TexCoords.y = (float)tex.mData[1];

								FbxVector4 norm;
								Mesh->GetPolygonVertexNormal(p, v, norm);
								NewVertex.Normal.x = (float)norm.mData[0];
								NewVertex.Normal.y = (float)norm.mData[1];
								NewVertex.Normal.z = (float)norm.mData[2];

								NewVertex.BoneIDBoneWeight = TempVertices[controlPointIndex].BoneIDBoneWeight;

								OutVertices.push_back(NewVertex);
								OutIndicies.back().Indices[v] = counter;

								counter++;
							}
						}

					

					}
					
					GetSubsetData(Mesh, OutIndicies, SubsetData);
					ProcessMaterials(pNode);
					this->AnimatedNodeMeshes.push_back(AnimatedMesh(NodeName, OutVertices, OutIndicies, SubsetData));
				}
			}
		}
		else
		{
			std::string nodeName = pNode->GetName();
			RemoveCharsFromString(nodeName, '.');
			if (pNode->GetParent())
			{
				this->ParentNodeID = pNode->GetParent()->GetUniqueID();
				this->NodeTransform = TransformNode(nodeName, parent);
			}
			else
			{
				this->ParentNodeID = -1;
				this->NodeTransform = TransformNode(nodeName, nullptr);
			}

			this->NodeID = pNode->GetUniqueID();
			this->NodeName = nodeName;
			this->NumberOfChilds = pNode->GetChildCount();
			this->NumberOfMeshes = 0;
		}
		for (int i = 0; i < pNode->GetChildCount(); i++)
		{
			FBXSceneNode child(pNode->GetChild(i)->GetName());
			child.ReadNodeMeshData(pNode->GetChild(i), &this->NodeTransform, pDepth + 1, skeleton, convertToLeftHanded, animations);
			this->ChildNodes.push_back(child);
		}
	}

	void SaveNodeDataToDataFile(DataFileContainer& MeshOutputFile, DataFileContainer& MaterialOutputFile,
	                            std::string MeshSection, std::string CollisionMeshOutputFilename,
	                            bool createCollisionMesh)
	{
		MeshOutputFile.AddIntValueByName(MeshSection, "Node_ID", NodeID);
		MeshOutputFile.AddIntValueByName(MeshSection, "Parent_Node_ID", ParentNodeID);
		MeshOutputFile.AddStringValueByName(MeshSection, "Node_Name", NodeName);
		//MeshOutputFile.AddIntValueByName(MeshSection, "Number_of_Childs", NumberOfChilds);
		MeshOutputFile.AddIntValueByName(MeshSection, "Number_of_Animated_Meshes", NumberOfMeshes);

		std::vector<float>pos;
		pos.push_back(DirectX::XMVectorGetX(NodeTransform.GetLocalPosition()));
		pos.push_back(DirectX::XMVectorGetY(NodeTransform.GetLocalPosition()));
		pos.push_back(DirectX::XMVectorGetZ(NodeTransform.GetLocalPosition()));

		std::vector<float>rotation;
		rotation.push_back(DirectX::XMVectorGetX(NodeTransform.GetLocalRotation()));
		rotation.push_back(DirectX::XMVectorGetY(NodeTransform.GetLocalRotation()));
		rotation.push_back(DirectX::XMVectorGetZ(NodeTransform.GetLocalRotation()));
		rotation.push_back(DirectX::XMVectorGetW(NodeTransform.GetLocalRotation()));

		std::vector<float>scale;
		scale.push_back(DirectX::XMVectorGetX(NodeTransform.GetLocalScale()));
		scale.push_back(DirectX::XMVectorGetY(NodeTransform.GetLocalScale()));
		scale.push_back(DirectX::XMVectorGetZ(NodeTransform.GetLocalScale()));


		MeshOutputFile.AddFloatVectorByName(MeshSection, "Node_Position", pos);
		MeshOutputFile.AddFloatVectorByName(MeshSection, "Node_Rotation", rotation);
		MeshOutputFile.AddFloatVectorByName(MeshSection, "Node_Scale", scale);

		for (int k = 0; k < NumberOfMeshes; k++)
		{
			std::vector<float> vertices;
			std::vector<int> indicies;
			std::vector<float> vertexPositions;

			if (createCollisionMesh)
			{
				MeshOutputFile.AddStringValueByName(MeshSection + ".AnimatedMesh_" + std::to_string(k), "Collision_Mesh_File", MeshSection + "_" + std::to_string(k) + "_" + CollisionMeshOutputFilename);
			}
			MeshOutputFile.AddStringValueByName(MeshSection + ".AnimatedMesh_" + std::to_string(k), "Mesh_Name", AnimatedNodeMeshes[k].MeshName);
			MeshOutputFile.AddIntValueByName(MeshSection + ".AnimatedMesh_" + std::to_string(k), "Subset_Count", 1);

			for (SubsetSurfaceData subset_data : AnimatedNodeMeshes[k].SubsetData)
			{
				MaterialOutputFile.AddBoolValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "render_forward", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].RenderForward);
				MaterialOutputFile.AddBoolValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "is_transparent", AnimatedNodeMeshes[k].SubsetData[0].SurfaceMaterial.is_transparent);
				MaterialOutputFile.AddBoolValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "use_alphachannel_transparency", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].useAlpaChannelTransparency);
				MaterialOutputFile.AddStringValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "albedo_map", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].AlbedoMapFilename);
				MaterialOutputFile.AddStringValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "roughness_map", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].RoughnessMapFilename);
				MaterialOutputFile.AddStringValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "metalness_map", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].MetalnessMapFilename);
				MaterialOutputFile.AddStringValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "normal_map", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].NormalMapFilename);
				MaterialOutputFile.AddFloatValueByName(mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName, "f0", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].f0);
			}


			for (int i = 0; i < AnimatedNodeMeshes[k].GetNumberOfVertices(); i++)
			{
				AnimatedVertex TempVertex = AnimatedNodeMeshes[k].GetVertexData(i);

				vertices.push_back(TempVertex.Position.x);
				vertices.push_back(TempVertex.Position.y);
				vertices.push_back(TempVertex.Position.z);

				vertices.push_back(TempVertex.TexCoords.x);
				vertices.push_back(TempVertex.TexCoords.y);

				vertices.push_back(TempVertex.Normal.x);
				vertices.push_back(TempVertex.Normal.y);
				vertices.push_back(TempVertex.Normal.z);

				for (auto bone_id_bone_weight : TempVertex.BoneIDBoneWeight)
				{
					if (bone_id_bone_weight.BoneId < 0)
					{
						vertices.push_back(-1);
						vertices.push_back(0);
					}
					else
					{
						vertices.push_back(bone_id_bone_weight.BoneId);
						vertices.push_back(bone_id_bone_weight.BoneWeight);
					}
				}

				vertexPositions.push_back(TempVertex.Position.x);
				vertexPositions.push_back(TempVertex.Position.y);
				vertexPositions.push_back(TempVertex.Position.z);
			}

			for (int i = 0; i < AnimatedNodeMeshes[k].GetNumberOfTriangleIndicies(); i++)
			{
				indicies.push_back(AnimatedNodeMeshes[k].GetTriangleIndexData(i).Indices[0]);
				indicies.push_back(AnimatedNodeMeshes[k].GetTriangleIndexData(i).Indices[1]);
				indicies.push_back(AnimatedNodeMeshes[k].GetTriangleIndexData(i).Indices[2]);
			}

			// Add the mesh vertices to the mesh output file.
			MeshOutputFile.AddFloatVectorByName(MeshSection + ".AnimatedMesh_" + std::to_string(k), "Vertices", vertices);
			MeshOutputFile.AddIntVectorByName(MeshSection + ".AnimatedMesh_" + std::to_string(k), "Indicies", indicies);

			int subsetIndex = 0;
			for (SubsetSurfaceData subset_data : AnimatedNodeMeshes[k].SubsetData)
			{
				std::string MeshSubsetSection = "_Subset_"+ std::to_string(subsetIndex);
				MeshOutputFile.AddIntValueByName(MeshSection + ".AnimatedMesh_" + std::to_string(k) + MeshSubsetSection, "Subset_Start", subset_data.Start);
				MeshOutputFile.AddIntValueByName(MeshSection + ".AnimatedMesh_" + std::to_string(k) + MeshSubsetSection, "Subset_Draw_Amount", subset_data.DrawAmount);
				MeshOutputFile.AddStringValueByName(MeshSection + ".AnimatedMesh_" + std::to_string(k) + MeshSubsetSection, "Subset_Material_Name", mMaterialLookUp[subset_data.SurfaceMaterial.materialIndex].matName);
			}

			if (createCollisionMesh)
			{
				std::cout << "Erstelle Collision Mesh " + MeshSection << std::endl;
				std::cout << "Dies kann einige Minuten in Anspruch nehmen.... " << std::endl;

				DecompResults* decompResults = Decomp(indicies, vertexPositions);

				// Erstelle convex decomposition mesh.
				btCompoundShape* mesh = decompResults->CompoundShape;

				// Speichere convex decomposition mesh
				btDefaultSerializer* serializer = new btDefaultSerializer();

				serializer->startSerialization();
				mesh->serializeSingleShape(serializer);
				serializer->finishSerialization();

				FILE* file;
				fopen_s(&file, (MeshSection + "_" + std::to_string(k) + "_" + CollisionMeshOutputFilename).c_str(), "wb");
				fwrite(serializer->getBufferPointer(), serializer->getCurrentBufferSize(), 1, file);
				fclose(file);

				// Free up used memory
				delete mesh;
				delete serializer;
			}

		}
		

		// Then do the same for each of its children
		for (int i = 0; i < NumberOfChilds; i++)
		{
			FBXSceneNode childNode = std::move(ChildNodes[i]);
			childNode.SaveNodeDataToDataFile(MeshOutputFile, MaterialOutputFile, childNode.NodeName, CollisionMeshOutputFilename, createCollisionMesh);
		}
	}
};

FBXSceneNode::FBXSceneNode(std::string nodeName)
	: NodeTransform(nodeName, nullptr),
	  NodeName(nodeName)
{
}

FBXSceneNode::FBXSceneNode(std::string nodeName, uint64_t parentNodeId, TransformNode* parentTransform, uint64_t nodeId,
                           DirectX::FXMVECTOR position, DirectX::FXMVECTOR quaternionRotation, DirectX::FXMVECTOR scale)
	: ParentNodeID(parentNodeId), NodeTransform(nodeName, parentTransform, position, quaternionRotation, scale),
	  NodeName(nodeName), NodeID(nodeId)
{
}

FBXSceneNode::FBXSceneNode(std::string nodeName, uint64_t parentNodeId, TransformNode* parentTransform, uint64_t nodeId)
	: ParentNodeID(parentNodeId), NodeTransform(nodeName, parentTransform),
	  NodeName(nodeName), NodeID(nodeId)
{
}

class FBXModel
{
public:
	FBXModel()
		: ModelRootNode("ModelRootNode")
	{
	}

	void ReadScene(FbxScene* pScene)
	{
		//FbxAxisSystem::DirectX.DeepConvertScene(pScene);
		bool shouldConvert = false;
		FbxAxisSystem SceneAxisSystem = pScene->GetGlobalSettings().GetAxisSystem();
		FbxAxisSystem OurAxisSystem(FbxAxisSystem::EPreDefinedAxisSystem::eDirectX);
		if (SceneAxisSystem != OurAxisSystem)
		{
			//shouldConvert = true;
			OurAxisSystem.ConvertScene(pScene);
			OurAxisSystem.DeepConvertScene(pScene);
		}

		FbxSystemUnit SceneSystemUnit = pScene->GetGlobalSettings().GetSystemUnit();
		if (pScene->GetGlobalSettings().GetSystemUnit() != FbxSystemUnit::m)
		{
			// Convert the scene to meters using the defined options.
			FbxSystemUnit::m.ConvertScene(pScene);
		}

		for (int animationStackIndex = 0; animationStackIndex < pScene->GetSrcObjectCount<FbxAnimStack>(); animationStackIndex++)
		{
			FBXAnimations.push_back(FBXAnimation());
			FbxAnimStack* currAnimStack = pScene->GetSrcObject<FbxAnimStack>(animationStackIndex);
			FbxString animStackName = currAnimStack->GetName();
			FBXAnimations.back().AnimationName = animStackName.Buffer();
			RemoveCharsFromString(FBXAnimations.back().AnimationName, '.');
			FbxTakeInfo* takeInfo = pScene->GetTakeInfo(animStackName);
			FbxTime start = takeInfo->mLocalTimeSpan.GetStart();
			FbxTime end = takeInfo->mLocalTimeSpan.GetStop();
			FBXAnimations.back().AnimationLength = end.GetMilliSeconds();
		}

		int i;
		FbxNode* lRootNode = pScene->GetRootNode();
		FBXSceneNode::ReadNodeSkeletonData(lRootNode, 0, -1, -1, Skeleton);
		ModelRootNode.ReadNodeMeshData(lRootNode, nullptr, 0, Skeleton, shouldConvert, FBXAnimations);

		std::map<std::string, AnimationEngineSide::BoneInfo> boneInfo;

		for (Joint m_joint : Skeleton.mJoints)
		{
			boneInfo[m_joint.mName].BoneId = m_joint.BoneId;
			boneInfo[m_joint.mName].BoneName = m_joint.mName;
			DirectX::XMStoreFloat4x4(&boneInfo[m_joint.mName].GlobalBindposeInverse, ConvertMatrix(m_joint.mGlobalBindposeInverse));
		}

		for (FBXAnimation animation : FBXAnimations)
		{
			
			std::vector<AnimationEngineSide::BoneAnimationKeyframes> BoneAnimationKeyframes;
			std::vector<FBXBoneAnimationKeyframes> AnimationKeyframes = animation.AnimationKeyframes;
			for (FBXBoneAnimationKeyframes animation_keyframe : AnimationKeyframes)
			{
				std::vector<AnimationEngineSide::KeyPosition> keyframesPosition;
				std::vector<AnimationEngineSide::KeyRotation> keyframesRotation;
				std::vector<AnimationEngineSide::KeyScale> keyframesScale;
				std::string boneName = animation_keyframe.BoneName;
				int boneId = animation_keyframe.BoneId;
				for (KeySRT bone_keyframe : animation_keyframe.boneKeyframes)
				{
					AnimationEngineSide::KeyPosition tempKeyPos;
					AnimationEngineSide::KeyRotation tempKeyRot;
					AnimationEngineSide::KeyScale tempKeyScale;
					ConvertMatrix(bone_keyframe.mGlobalTransform, tempKeyScale.scale, tempKeyRot.orientation, tempKeyPos.position);
					tempKeyPos.timeStamp = bone_keyframe.Timestamp;
					tempKeyRot.timeStamp = bone_keyframe.Timestamp;
					tempKeyScale.timeStamp = bone_keyframe.Timestamp;
					keyframesPosition.push_back(tempKeyPos);
					keyframesRotation.push_back(tempKeyRot);
					keyframesScale.push_back(tempKeyScale);
				}

				BoneAnimationKeyframes.push_back(AnimationEngineSide::BoneAnimationKeyframes(boneName, boneId, keyframesPosition, keyframesRotation, keyframesScale));
			}

			AnimationsEngineSide.push_back(AnimationEngineSide::Animation(animation.AnimationName, animation.AnimationLength, 24, BoneAnimationKeyframes, boneInfo, &ModelRootNode.NodeTransform));
		}
		
	}

	void SaveScene(std::string modelName, bool createCollisionMeshes)
	{
		DataFileContainer MeshOutputFile;
		DataFileContainer MaterialOutputFile;
		DataFileContainer AnimationOutputFile;
		MaterialOutputFile.ClearDataContainer();
		MeshOutputFile.ClearDataContainer();

		std::string MeshOutputFilename = modelName + ".smo";
		std::string MaterialOutputFilename = modelName + ".matlib";
		std::string CollisionMeshOutputFilename = modelName + ".bcs";
		std::string AnimationOutputFilename = modelName + ".anim";

		MeshOutputFile.AddStringValueByName("Model_Information", "MaterialLibrary", MaterialOutputFilename);

		ModelRootNode.SaveNodeDataToDataFile(MeshOutputFile, MaterialOutputFile, ModelRootNode.NodeName, CollisionMeshOutputFilename, createCollisionMeshes);
		SaveAnimations(AnimationOutputFile);

		MeshOutputFile.SaveDataContainerToBinaryFile(MeshOutputFilename);
		MaterialOutputFile.SaveDataContainerToFile(MaterialOutputFilename, "Material Library");
		AnimationOutputFile.SaveDataContainerToBinaryFile(AnimationOutputFilename);

		MeshOutputFile.SaveDataContainerToFile(MeshOutputFilename + ".test", "MeshFileTest");
		AnimationOutputFile.SaveDataContainerToFile(AnimationOutputFilename + ".test", "AnimationFileTest");
	}

private:

	void SaveAnimations(DataFileContainer& AnimationOutputFile)
	{
		AnimationOutputFile.AddIntValueByName("Animations", "NumberOfAnimations", AnimationsEngineSide.size());
		for (auto animation : AnimationsEngineSide)
		{
			AnimationOutputFile.AddFloatValueByName("Animations." + animation.GetAnimationName(), "AnimationDuration", animation.GetDuration());
			AnimationOutputFile.AddFloatValueByName("Animations." + animation.GetAnimationName(), "AnimationTicksPerSecond", animation.GetTicksPerSecond());
			std::vector<AnimationEngineSide::BoneAnimationKeyframes> boneKeyframes = animation.GetBonesAnimationInformation();
			for (AnimationEngineSide::BoneAnimationKeyframes bone_keyframe : boneKeyframes)
			{
				std::vector<AnimationEngineSide::KeyPosition> bonePositionKeyframes = bone_keyframe.GetPositionKeyframes();
				int count = 0;
				for (AnimationEngineSide::KeyPosition bone_position_keyframe : bonePositionKeyframes)
				{
					AnimationOutputFile.AddFloatValueByName("Animations." + animation.GetAnimationName() + "." + bone_keyframe.GetBoneName() + ".PositionKeys.PositionKey_" + std::to_string(count), "TimeStamp", bone_position_keyframe.timeStamp);
					std::vector<float>positionKey;
					positionKey.push_back(bone_position_keyframe.position.x);
					positionKey.push_back(bone_position_keyframe.position.y);
					positionKey.push_back(bone_position_keyframe.position.z);
					AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + bone_keyframe.GetBoneName() + ".PositionKeys.PositionKey_" + std::to_string(count), "Position", positionKey);
					count++;
				}

				std::vector<AnimationEngineSide::KeyRotation> boneRotationKeyframes = bone_keyframe.GetRotationKeyframes();
				count = 0;
				for (AnimationEngineSide::KeyRotation bone_position_keyframe : boneRotationKeyframes)
				{
					AnimationOutputFile.AddFloatValueByName("Animations." + animation.GetAnimationName() + "." + bone_keyframe.GetBoneName() + ".RotationKeys.RotationKey_" + std::to_string(count), "TimeStamp", bone_position_keyframe.timeStamp);
					std::vector<float>positionKey;
					positionKey.push_back(bone_position_keyframe.orientation.x);
					positionKey.push_back(bone_position_keyframe.orientation.y);
					positionKey.push_back(bone_position_keyframe.orientation.z);
					positionKey.push_back(bone_position_keyframe.orientation.w);
					AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + bone_keyframe.GetBoneName() + ".RotationKeys.RotationKey_" + std::to_string(count), "Rotation", positionKey);
					count++;
				}

				std::vector<AnimationEngineSide::KeyScale> boneScaleKeyframes = bone_keyframe.GetScaleKeyframes();
				count = 0;
				for (AnimationEngineSide::KeyScale bone_position_keyframe : boneScaleKeyframes)
				{
					AnimationOutputFile.AddFloatValueByName("Animations." + animation.GetAnimationName() + "." + bone_keyframe.GetBoneName() + ".ScaleKeys.ScaleKey_" + std::to_string(count), "TimeStamp", bone_position_keyframe.timeStamp);
					std::vector<float>positionKey;
					positionKey.push_back(bone_position_keyframe.scale.x);
					positionKey.push_back(bone_position_keyframe.scale.y);
					positionKey.push_back(bone_position_keyframe.scale.z);
					AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + bone_keyframe.GetBoneName() + ".ScaleKeys.ScaleKey_" + std::to_string(count), "Scale", positionKey);
					count++;
				}


			}


			std::map<std::string, AnimationEngineSide::BoneInfo> boneInfoMap = animation.GetBoneIDMap();

			int count = 0;
			for (auto boneInfo : boneInfoMap)
			{
				if (boneInfo.second.BoneId < 0)
				{
					ASSERT(0);
				}
				AnimationOutputFile.AddIntValueByName("Animations." + animation.GetAnimationName() + "." + boneInfo.first, "BoneId", boneInfo.second.BoneId);
				std::vector<float>positionKey;
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._11);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._12);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._13);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._14);
				AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + boneInfo.first + ".BoneInfoMap", "BoneOffsetMatrixRowOne", positionKey);
				positionKey.clear();
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._21);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._22);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._23);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._24);
				AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + boneInfo.first + ".BoneInfoMap", "BoneOffsetMatrixRowTwo", positionKey);
				positionKey.clear();
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._31);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._32);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._33);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._34);
				AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + boneInfo.first + ".BoneInfoMap", "BoneOffsetMatrixRowThree", positionKey);
				positionKey.clear();
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._41);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._42);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._43);
				positionKey.push_back(boneInfo.second.GlobalBindposeInverse._44);
				AnimationOutputFile.AddFloatVectorByName("Animations." + animation.GetAnimationName() + "." + boneInfo.first + ".BoneInfoMap", "BoneOffsetMatrixRowFour", positionKey);
				count++;
			}
		}
	}

	std::vector<FBXAnimation>FBXAnimations;
	std::vector<AnimationEngineSide::Animation>AnimationsEngineSide;
	FBXSceneNode ModelRootNode;
	Skeleton Skeleton;
};

#ifdef IOS_REF
#undef  IOS_REF
#define IOS_REF (*(pManager->GetIOSettings()))
#endif

void InitializeSdkObjects(FbxManager*& pManager, FbxScene*& pScene)
{
	//The first thing to do is to create the FBX Manager which is the object allocator for almost all the classes in the SDK
	pManager = FbxManager::Create();
	if (!pManager)
	{
		FBXSDK_printf("Error: Unable to create FBX Manager!\n");
		exit(1);
	}
	else
		FBXSDK_printf("Autodesk FBX SDK version %s\n", pManager->GetVersion());

	//Create an IOSettings object. This object holds all import/export settings.
	FbxIOSettings* ios = FbxIOSettings::Create(pManager, IOSROOT);
	pManager->SetIOSettings(ios);

	//Load plugins from the executable directory (optional)
	FbxString lPath = FbxGetApplicationDirectory();
	pManager->LoadPluginsDirectory(lPath.Buffer());

	//Create an FBX scene. This object holds most objects imported/exported from/to files.
	pScene = FbxScene::Create(pManager, "My Scene");
	if (!pScene)
	{
		FBXSDK_printf("Error: Unable to create FBX scene!\n");
		exit(1);
	}
}

void DestroySdkObjects(FbxManager* pManager, bool pExitStatus)
{
	//Delete the FBX Manager. All the objects that have been allocated using the FBX Manager and that haven't been explicitly destroyed are also automatically destroyed.
	if (pManager) pManager->Destroy();
	if (pExitStatus)
		FBXSDK_printf("Program Success!\n");
}

bool LoadScene(FbxManager* pManager, FbxDocument* pScene, const char* pFilename)
{
	int lFileMajor, lFileMinor, lFileRevision;
	int lSDKMajor, lSDKMinor, lSDKRevision;
	//int lFileFormat = -1;
	int lAnimStackCount;
	bool lStatus;
	char lPassword[1024];

	// Get the file version number generate by the FBX SDK.
	FbxManager::GetFileFormatVersion(lSDKMajor, lSDKMinor, lSDKRevision);

	// Create an importer.
	FbxImporter* lImporter = FbxImporter::Create(pManager, "");

	// Initialize the importer by providing a filename.
	const bool lImportStatus = lImporter->Initialize(pFilename, -1, pManager->GetIOSettings());
	lImporter->GetFileVersion(lFileMajor, lFileMinor, lFileRevision);

	if (!lImportStatus)
	{
		FbxString error = lImporter->GetStatus().GetErrorString();
		FBXSDK_printf("Call to FbxImporter::Initialize() failed.\n");
		FBXSDK_printf("Error returned: %s\n\n", error.Buffer());

		if (lImporter->GetStatus().GetCode() == FbxStatus::eInvalidFileVersion)
		{
			FBXSDK_printf("FBX file format version for this FBX SDK is %d.%d.%d\n", lSDKMajor, lSDKMinor, lSDKRevision);
			FBXSDK_printf("FBX file format version for file '%s' is %d.%d.%d\n\n", pFilename, lFileMajor, lFileMinor,
			              lFileRevision);
		}

		return false;
	}

	FBXSDK_printf("FBX file format version for this FBX SDK is %d.%d.%d\n", lSDKMajor, lSDKMinor, lSDKRevision);

	if (lImporter->IsFBX())
	{
		FBXSDK_printf("FBX file format version for file '%s' is %d.%d.%d\n\n", pFilename, lFileMajor, lFileMinor,
		              lFileRevision);

		// From this point, it is possible to access animation stack information without
		// the expense of loading the entire file.

		FBXSDK_printf("Animation Stack Information\n");

		lAnimStackCount = lImporter->GetAnimStackCount();

		FBXSDK_printf("    Number of Animation Stacks: %d\n", lAnimStackCount);
		FBXSDK_printf("    Current Animation Stack: \"%s\"\n", lImporter->GetActiveAnimStackName().Buffer());
		FBXSDK_printf("\n");

		for (int i = 0; i < lAnimStackCount; i++)
		{
			FbxTakeInfo* lTakeInfo = lImporter->GetTakeInfo(i);

			FBXSDK_printf("    Animation Stack %d\n", i);
			FBXSDK_printf("         Name: \"%s\"\n", lTakeInfo->mName.Buffer());
			FBXSDK_printf("         Description: \"%s\"\n", lTakeInfo->mDescription.Buffer());

			// Change the value of the import name if the animation stack should be imported 
			// under a different name.
			FBXSDK_printf("         Import Name: \"%s\"\n", lTakeInfo->mImportName.Buffer());

			// Set the value of the import state to false if the animation stack should be not
			// be imported. 
			FBXSDK_printf("         Import State: %s\n", lTakeInfo->mSelect ? "true" : "false");
			FBXSDK_printf("\n");
		}

		// Set the import states. By default, the import states are always set to 
		// true. The code below shows how to change these states.
		IOS_REF.SetBoolProp(IMP_FBX_MATERIAL, true);
		IOS_REF.SetBoolProp(IMP_FBX_TEXTURE, true);
		IOS_REF.SetBoolProp(IMP_FBX_LINK, true);
		IOS_REF.SetBoolProp(IMP_FBX_SHAPE, true);
		IOS_REF.SetBoolProp(IMP_FBX_GOBO, true);
		IOS_REF.SetBoolProp(IMP_FBX_ANIMATION, true);
		IOS_REF.SetBoolProp(IMP_FBX_GLOBAL_SETTINGS, true);
	}

	// Import the scene.
	lStatus = lImporter->Import(pScene);
	if (lStatus == false && lImporter->GetStatus() == FbxStatus::ePasswordError)
	{
		FBXSDK_printf("Please enter password: ");

		lPassword[0] = '\0';

		FBXSDK_CRT_SECURE_NO_WARNING_BEGIN
		scanf("%s", lPassword);
		FBXSDK_CRT_SECURE_NO_WARNING_END

		FbxString lString(lPassword);

		IOS_REF.SetStringProp(IMP_FBX_PASSWORD, lString);
		IOS_REF.SetBoolProp(IMP_FBX_PASSWORD_ENABLE, true);

		lStatus = lImporter->Import(pScene);

		if (lStatus == false && lImporter->GetStatus() == FbxStatus::ePasswordError)
		{
			FBXSDK_printf("\nPassword is wrong, import aborted.\n");
		}
	}

	if (!lStatus || (lImporter->GetStatus() != FbxStatus::eSuccess))
	{
		FBXSDK_printf("********************************************************************************\n");
		if (lStatus)
		{
			FBXSDK_printf("WARNING:\n");
			FBXSDK_printf("   The importer was able to read the file but with errors.\n");
			FBXSDK_printf("   Loaded scene may be incomplete.\n\n");
		}
		else
		{
			FBXSDK_printf("Importer failed to load the file!\n\n");
		}

		if (lImporter->GetStatus() != FbxStatus::eSuccess)
			FBXSDK_printf("   Last error message: %s\n", lImporter->GetStatus().GetErrorString());

		FbxArray<FbxString*> history;
		lImporter->GetStatus().GetErrorStringHistory(history);
		if (history.GetCount() > 1)
		{
			FBXSDK_printf("   Error history stack:\n");
			for (int i = 0; i < history.GetCount(); i++)
			{
				FBXSDK_printf("      %s\n", history[i]->Buffer());
			}
		}
		FbxArrayDelete<FbxString*>(history);
		FBXSDK_printf("********************************************************************************\n");
	}

	FbxGeometryConverter clsConverter(pManager);
	clsConverter.Triangulate(pScene->GetScene(), true);

	// Destroy the importer.
	lImporter->Destroy();

	return lStatus;
}





void main()
{
	bool result;
	std::string MeshFilename;
	FbxManager* lSdkManager = NULL;
	FbxScene* lScene = NULL;
	bool lResult;

	InitializeSdkObjects(lSdkManager, lScene);

	std::cout << "Geben Sie den Dateinamen des umzuwandelnen Model ein: ";
	MeshFilename = "HipHopDancing.fbx";

	FbxString lFilePath("");
	lFilePath = MeshFilename.c_str();
	FBXModel fbxModel;
	if (lFilePath.IsEmpty())
	{
		lResult = false;
		FBXSDK_printf("\n\nUsage: ImportScene <FBX file name>\n\n");
	}
	else
	{
		FBXSDK_printf("\n\nFile: %s\n\n", lFilePath.Buffer());
		lResult = LoadScene(lSdkManager, lScene, lFilePath.Buffer());
	}
	if (lResult == false)
	{
		FBXSDK_printf("\n\nAn error occurred while loading the scene...");
	}
	else
	{
		fbxModel.ReadScene(lScene);

		std::filesystem::path p(MeshFilename);
		fbxModel.SaveScene(p.stem().string(), false);
	}
	DestroySdkObjects(lSdkManager, false);

	system("pause");
}
