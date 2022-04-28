#pragma once

#include <fbxsdk.h>

#include "Mesh.h"

#include "AnimationStuff.h"



inline DirectX::XMMATRIX ConvertMatrix(FbxAMatrix inMat)
{
	FbxVector4 scaling = inMat.GetS();
	FbxQuaternion rotation = inMat.GetQ();
	FbxVector4 translation = inMat.GetT();

	DirectX::XMFLOAT3 Scaling;
	DirectX::XMFLOAT4 Rotation;
	DirectX::XMFLOAT3 Translation;

	Scaling.x = scaling.mData[0];
	Scaling.y = scaling.mData[1];
	Scaling.z = scaling.mData[2];

	Rotation.x = rotation.mData[0];
	Rotation.y = rotation.mData[1];
	Rotation.z = rotation.mData[2];
	Rotation.w = rotation.mData[3];

	Translation.x = translation.mData[0];
	Translation.y = translation.mData[1];
	Translation.z = translation.mData[2];

	const DirectX::XMMATRIX ScaleMatrix = DirectX::XMMatrixScaling(Scaling.x, Scaling.y, Scaling.z);

	const DirectX::XMMATRIX RotationMatrix = DirectX::XMMatrixRotationQuaternion(DirectX::XMLoadFloat4(&Rotation));

	const DirectX::XMMATRIX TranslateMatrix = DirectX::XMMatrixTranslation(Translation.x, Translation.y, Translation.z);


	DirectX::XMMATRIX srt = ScaleMatrix * RotationMatrix * TranslateMatrix;

	return srt;
}

inline void fbxamatrix_to_xmfloat4x4(FbxAMatrix& fbxamatrix, DirectX::XMFLOAT4X4& xmfloat4x4)
{
	//fbxamatrix = fbxamatrix.Transpose();
	xmfloat4x4 = DirectX::XMFLOAT4X4((float)fbxamatrix.GetRow(0)[0],
													(float)fbxamatrix.GetRow(0)[1], (float)fbxamatrix.GetRow(0)[2],
													(float)fbxamatrix.GetRow(0)[3], (float)fbxamatrix.GetRow(1)[0],
													(float)fbxamatrix.GetRow(1)[1], (float)fbxamatrix.GetRow(1)[2],
													(float)fbxamatrix.GetRow(1)[3], (float)fbxamatrix.GetRow(2)[0],
													(float)fbxamatrix.GetRow(2)[1], (float)fbxamatrix.GetRow(2)[2],
													(float)fbxamatrix.GetRow(2)[3], (float)fbxamatrix.GetRow(3)[0],
													(float)fbxamatrix.GetRow(3)[1], (float)fbxamatrix.GetRow(3)[2],
													(float)fbxamatrix.GetRow(3)[3]);;
}

inline void ConvertMatrix(FbxAMatrix inMat, DirectX::XMFLOAT3& OutScaling, DirectX::XMFLOAT4& OutRotation, DirectX::XMFLOAT3& OutTranslation)
{
	FbxVector4 scaling = inMat.GetS();
	FbxQuaternion rotation = inMat.GetQ();
	FbxVector4 translation = inMat.GetT();

	OutScaling.x = scaling.mData[0];
	OutScaling.y = scaling.mData[1];
	OutScaling.z = scaling.mData[2];

	OutRotation.x = rotation.mData[0];
	OutRotation.y = rotation.mData[1];
	OutRotation.z = rotation.mData[2];
	OutRotation.w = rotation.mData[3];

	OutTranslation.x = translation.mData[0];
	OutTranslation.y = translation.mData[1];
	OutTranslation.z = translation.mData[2];
}

struct BoneInfoNew
{
	BoneInfoNew() = default;

	BoneInfoNew(const BoneInfoNew& other)
		: BoneId(other.BoneId),
		  BoneOffsetMatrix(other.BoneOffsetMatrix)
	{
	}

	BoneInfoNew(BoneInfoNew&& other) noexcept
		: BoneId(other.BoneId),
		  BoneOffsetMatrix(std::move(other.BoneOffsetMatrix))
	{
	}

	BoneInfoNew& operator=(const BoneInfoNew& other)
	{
		if (this == &other)
			return *this;
		BoneId = other.BoneId;
		BoneOffsetMatrix = other.BoneOffsetMatrix;
		return *this;
	}

	BoneInfoNew& operator=(BoneInfoNew&& other) noexcept
	{
		if (this == &other)
			return *this;
		BoneId = other.BoneId;
		BoneOffsetMatrix = std::move(other.BoneOffsetMatrix);
		return *this;
	}

	/*id is index in finalBoneMatrices*/
	int BoneId;

	/*offset matrix transforms vertex from model space to bone space*/
	DirectX::XMFLOAT4X4 BoneOffsetMatrix;

};

struct BoneIdBoneWeight
{
	int BoneId;
	float BoneWeight;
};

struct AnimatedVertex
{
	AnimatedVertex(const AnimatedVertex& other)
		: Position(other.Position),
		  Normal(other.Normal),
		  TexCoords(other.TexCoords),
		  BoneIDBoneWeight(other.BoneIDBoneWeight)
	{
	}

	AnimatedVertex(AnimatedVertex&& other) noexcept
		: Position(std::move(other.Position)),
		  Normal(std::move(other.Normal)),
		  TexCoords(std::move(other.TexCoords)),
		  BoneIDBoneWeight(std::move(other.BoneIDBoneWeight))
	{
	}

	AnimatedVertex& operator=(const AnimatedVertex& other)
	{
		if (this == &other)
			return *this;
		Position = other.Position;
		Normal = other.Normal;
		TexCoords = other.TexCoords;
		BoneIDBoneWeight = other.BoneIDBoneWeight;
		return *this;
	}

	AnimatedVertex& operator=(AnimatedVertex&& other) noexcept
	{
		if (this == &other)
			return *this;
		Position = std::move(other.Position);
		Normal = std::move(other.Normal);
		TexCoords = std::move(other.TexCoords);
		BoneIDBoneWeight = std::move(other.BoneIDBoneWeight);
		return *this;
	}

	AnimatedVertex();

	AnimatedVertex(DirectX::XMVECTOR VertexPosition, DirectX::XMVECTOR VertexNormal, DirectX::XMVECTOR VertexTexCoords, std::list<BoneIdBoneWeight> BoneIDBoneWeight);

	DirectX::XMFLOAT3 Position;
	DirectX::XMFLOAT3 Normal;
	DirectX::XMFLOAT2 TexCoords;
	std::list<BoneIdBoneWeight> BoneIDBoneWeight;
};

struct AnimatedMesh
{
	AnimatedMesh(const AnimatedMesh& other)
		: MeshName(other.MeshName),
		  Vertices(other.Vertices),
		  Indices(other.Indices),
		  SubsetData(other.SubsetData)
	{
	}

	AnimatedMesh(AnimatedMesh&& other) noexcept
		: MeshName(std::move(other.MeshName)),
		  Vertices(std::move(other.Vertices)),
		  Indices(std::move(other.Indices)),
		  SubsetData(std::move(other.SubsetData))
	{
	}

	AnimatedMesh& operator=(const AnimatedMesh& other)
	{
		if (this == &other)
			return *this;
		MeshName = other.MeshName;
		Vertices = other.Vertices;
		Indices = other.Indices;
		SubsetData = other.SubsetData;
		return *this;
	}

	AnimatedMesh& operator=(AnimatedMesh&& other) noexcept
	{
		if (this == &other)
			return *this;
		MeshName = std::move(other.MeshName);
		Vertices = std::move(other.Vertices);
		Indices = std::move(other.Indices);
		SubsetData = std::move(other.SubsetData);
		return *this;
	}


	AnimatedMesh(std::vector<AnimatedVertex> vertices, std::vector<TriangleVertexIndicies> indices, std::vector<SubsetSurfaceData> subsetData);

	AnimatedMesh(std::string meshName, std::vector<AnimatedVertex> vertices,
		std::vector<TriangleVertexIndicies> indices, std::vector<SubsetSurfaceData> subsetData);

	AnimatedVertex GetVertexData(int VertexIndex);
	int GetNumberOfVertices();

	TriangleVertexIndicies GetTriangleIndexData(int TriangleIndex);
	int GetNumberOfTriangleIndicies();

	std::string MeshName;

	std::vector<AnimatedVertex> Vertices;
	std::vector<TriangleVertexIndicies> Indices;
	std::vector<SubsetSurfaceData> SubsetData;

	


	//static void SetVertexBoneDataToDefault(AnimatedVertex& vertex);
};

namespace AnimationEngineSide
{
	struct BoneInfo
	{
		BoneInfo() = default;
		BoneInfo(const BoneInfo& other)
			: BoneId(other.BoneId),
			  BoneName(other.BoneName),
			  GlobalBindposeInverse(other.GlobalBindposeInverse)
		{
		}
		BoneInfo(BoneInfo&& other) noexcept
			: BoneId(other.BoneId),
			  BoneName(std::move(other.BoneName)),
			  GlobalBindposeInverse(std::move(other.GlobalBindposeInverse))
		{
		}
		BoneInfo& operator=(const BoneInfo& other)
		{
			if (this == &other)
				return *this;
			BoneId = other.BoneId;
			BoneName = other.BoneName;
			GlobalBindposeInverse = other.GlobalBindposeInverse;
			return *this;
		}
		BoneInfo& operator=(BoneInfo&& other) noexcept
		{
			if (this == &other)
				return *this;
			BoneId = other.BoneId;
			BoneName = std::move(other.BoneName);
			GlobalBindposeInverse = std::move(other.GlobalBindposeInverse);
			return *this;
		}
		/*id is index in finalBoneMatrices*/
		int BoneId;

		std::string BoneName;

		/*offset matrix transforms vertex from model space to bone space*/
		DirectX::XMFLOAT4X4 GlobalBindposeInverse;

	};


	struct KeyPosition
	{
		KeyPosition(const KeyPosition& other)
			: position(other.position),
			  timeStamp(other.timeStamp)
		{
		}
		KeyPosition(KeyPosition&& other) noexcept
			: position(std::move(other.position)),
			  timeStamp(other.timeStamp)
		{
		}
		KeyPosition& operator=(const KeyPosition& other)
		{
			if (this == &other)
				return *this;
			position = other.position;
			timeStamp = other.timeStamp;
			return *this;
		}
		KeyPosition& operator=(KeyPosition&& other) noexcept
		{
			if (this == &other)
				return *this;
			position = std::move(other.position);
			timeStamp = other.timeStamp;
			return *this;
		}
		KeyPosition() = default;
		DirectX::XMFLOAT3 position;
		float timeStamp;
	};

	struct KeyRotation
	{
		KeyRotation() = default;
		KeyRotation(const KeyRotation& other)
			: orientation(other.orientation),
			  timeStamp(other.timeStamp)
		{
		}
		KeyRotation(KeyRotation&& other) noexcept
			: orientation(std::move(other.orientation)),
			  timeStamp(other.timeStamp)
		{
		}
		KeyRotation& operator=(const KeyRotation& other)
		{
			if (this == &other)
				return *this;
			orientation = other.orientation;
			timeStamp = other.timeStamp;
			return *this;
		}
		KeyRotation& operator=(KeyRotation&& other) noexcept
		{
			if (this == &other)
				return *this;
			orientation = std::move(other.orientation);
			timeStamp = other.timeStamp;
			return *this;
		}
		DirectX::XMFLOAT4 orientation;
		float timeStamp;
	};

	struct KeyScale
	{
		KeyScale(const KeyScale& other)
			: scale(other.scale),
			  timeStamp(other.timeStamp)
		{
		}
		KeyScale(KeyScale&& other) noexcept
			: scale(std::move(other.scale)),
			  timeStamp(other.timeStamp)
		{
		}
		KeyScale& operator=(const KeyScale& other)
		{
			if (this == &other)
				return *this;
			scale = other.scale;
			timeStamp = other.timeStamp;
			return *this;
		}
		KeyScale& operator=(KeyScale&& other) noexcept
		{
			if (this == &other)
				return *this;
			scale = std::move(other.scale);
			timeStamp = other.timeStamp;
			return *this;
		}
		KeyScale() = default;
		DirectX::XMFLOAT3 scale;
		float timeStamp;
	};

	class BoneAnimationKeyframes
	{
	public:
		BoneAnimationKeyframes(const BoneAnimationKeyframes& other)
			: m_Positions(other.m_Positions),
			  m_Rotations(other.m_Rotations),
			  m_Scales(other.m_Scales),
			  m_NumPositions(other.m_NumPositions),
			  m_NumRotations(other.m_NumRotations),
			  m_NumScalings(other.m_NumScalings),
			  m_LocalTransform(other.m_LocalTransform),
			  m_Name(other.m_Name),
			  m_ID(other.m_ID)
		{
		}
		BoneAnimationKeyframes(BoneAnimationKeyframes&& other) noexcept
			: m_Positions(std::move(other.m_Positions)),
			  m_Rotations(std::move(other.m_Rotations)),
			  m_Scales(std::move(other.m_Scales)),
			  m_NumPositions(other.m_NumPositions),
			  m_NumRotations(other.m_NumRotations),
			  m_NumScalings(other.m_NumScalings),
			  m_LocalTransform(std::move(other.m_LocalTransform)),
			  m_Name(std::move(other.m_Name)),
			  m_ID(other.m_ID)
		{
		}
		BoneAnimationKeyframes& operator=(const BoneAnimationKeyframes& other)
		{
			if (this == &other)
				return *this;
			m_Positions = other.m_Positions;
			m_Rotations = other.m_Rotations;
			m_Scales = other.m_Scales;
			m_NumPositions = other.m_NumPositions;
			m_NumRotations = other.m_NumRotations;
			m_NumScalings = other.m_NumScalings;
			m_LocalTransform = other.m_LocalTransform;
			m_Name = other.m_Name;
			m_ID = other.m_ID;
			return *this;
		}
		BoneAnimationKeyframes& operator=(BoneAnimationKeyframes&& other) noexcept
		{
			if (this == &other)
				return *this;
			m_Positions = std::move(other.m_Positions);
			m_Rotations = std::move(other.m_Rotations);
			m_Scales = std::move(other.m_Scales);
			m_NumPositions = other.m_NumPositions;
			m_NumRotations = other.m_NumRotations;
			m_NumScalings = other.m_NumScalings;
			m_LocalTransform = std::move(other.m_LocalTransform);
			m_Name = std::move(other.m_Name);
			m_ID = other.m_ID;
			return *this;
		}
		BoneAnimationKeyframes(const std::string& name, int ID, std::vector<KeyPosition> keyPositions, std::vector<KeyRotation> keyRotations, std::vector<KeyScale> keyScales)
			:
			m_Name(name),
			m_ID(ID)
		{
			DirectX::XMStoreFloat4x4(&m_LocalTransform, DirectX::XMMatrixIdentity());

			m_NumPositions = keyPositions.size();
			m_Positions = keyPositions;

			m_NumRotations = keyRotations.size();
			m_Rotations = keyRotations;

			m_NumScalings = keyScales.size();
			m_Scales = keyScales;
		}

		void Update(float animationTime)
		{
			DirectX::XMMATRIX translation = InterpolatePosition(animationTime);
			DirectX::XMMATRIX rotation = InterpolateRotation(animationTime);
			DirectX::XMMATRIX scale = InterpolateScaling(animationTime);
			DirectX::XMStoreFloat4x4(&m_LocalTransform, scale * rotation * translation);
		}
		DirectX::XMMATRIX GetLocalTransform() { return DirectX::XMLoadFloat4x4(&m_LocalTransform); }
		std::string GetBoneName() const { return m_Name; }
		int GetBoneID() { return m_ID; }

		std::vector<KeyPosition> GetPositionKeyframes()
		{
			return m_Positions;
		}

		std::vector<KeyRotation> GetRotationKeyframes()
		{
			return m_Rotations;
		}

		std::vector<KeyScale> GetScaleKeyframes()
		{
			return m_Scales;
		}


		int GetPositionIndex(float animationTime)
		{
			for (int index = 0; index < m_NumPositions - 1; ++index)
			{
				if (animationTime < m_Positions[index + 1].timeStamp)
					return index;
			}
			return -1;
		}

		int GetRotationIndex(float animationTime)
		{
			for (int index = 0; index < m_NumRotations - 1; ++index)
			{
				if (animationTime < m_Rotations[index + 1].timeStamp)
					return index;
			}
			return -1;
		}

		int GetScaleIndex(float animationTime)
		{
			for (int index = 0; index < m_NumScalings - 1; ++index)
			{
				if (animationTime < m_Scales[index + 1].timeStamp)
					return index;
			}
			return -1;
		}


	private:

		float GetScaleFactor(float lastTimeStamp, float nextTimeStamp, float animationTime)
		{
			float scaleFactor = 0.0f;
			float midWayLength = animationTime - lastTimeStamp;
			float framesDiff = nextTimeStamp - lastTimeStamp;
			scaleFactor = midWayLength / framesDiff;
			return scaleFactor;
		}

		DirectX::XMMATRIX InterpolatePosition(float animationTime)
		{
			if (1 == m_NumPositions)
				return DirectX::XMMatrixTranslationFromVector(DirectX::XMLoadFloat3(&m_Positions[0].position));

			int p0Index = GetPositionIndex(animationTime);
			int p1Index = p0Index + 1;
			if (p0Index == -1)
			{
				return DirectX::XMMatrixTranslation(0.0f, 0.0f, 0.0f);
			}
			else
			{
				float scaleFactor = GetScaleFactor(m_Positions[p0Index].timeStamp,
				                                   m_Positions[p1Index].timeStamp, animationTime);
				const DirectX::XMVECTOR finalPosition = DirectX::XMVectorLerp(DirectX::XMLoadFloat3(&m_Positions[p0Index].position), DirectX::XMLoadFloat3(&m_Positions[p1Index].position)
				                                                              , scaleFactor);
				return  DirectX::XMMatrixTranslationFromVector(finalPosition);
			}
		}

		DirectX::XMMATRIX InterpolateRotation(float animationTime)
		{
			if (1 == m_NumRotations)
			{
				auto rotation = DirectX::XMVector4Normalize(DirectX::XMLoadFloat4(&m_Rotations[0].orientation));
				return DirectX::XMMatrixRotationQuaternion(rotation);
			}

			int p0Index = GetRotationIndex(animationTime);
			int p1Index = p0Index + 1;

			if (p0Index == -1)
			{
				return DirectX::XMMatrixRotationQuaternion(DirectX::XMQuaternionIdentity());
			}
			else
			{
				float scaleFactor = GetScaleFactor(m_Rotations[p0Index].timeStamp,
				                                   m_Rotations[p1Index].timeStamp, animationTime);
				DirectX::XMVECTOR finalRotation = DirectX::XMQuaternionSlerp(DirectX::XMLoadFloat4(&m_Rotations[p0Index].orientation), DirectX::XMLoadFloat4(&m_Rotations[p1Index].orientation)
				                                                             , scaleFactor);
				finalRotation = DirectX::XMQuaternionNormalize(finalRotation);
				return DirectX::XMMatrixRotationQuaternion(finalRotation);
			}

		}

		DirectX::XMMATRIX InterpolateScaling(float animationTime)
		{
			if (1 == m_NumScalings)
				return DirectX::XMMatrixScalingFromVector(DirectX::XMLoadFloat3(&m_Scales[0].scale));

			int p0Index = GetScaleIndex(animationTime);
			int p1Index = p0Index + 1;


			if (p0Index == -1)
			{
				return DirectX::XMMatrixScaling(1.0f, 1.0f, 1.0f);
			}
			else
			{
				float scaleFactor = GetScaleFactor(m_Scales[p0Index].timeStamp,
				                                   m_Scales[p1Index].timeStamp, animationTime);
				const DirectX::XMVECTOR finalScale = DirectX::XMVectorLerp(DirectX::XMLoadFloat3(&m_Scales[p0Index].scale), DirectX::XMLoadFloat3(&m_Scales[p1Index].scale)
				                                                           , scaleFactor);
				return  DirectX::XMMatrixScalingFromVector(finalScale);
			}

		}

		std::vector<KeyPosition> m_Positions;
		std::vector<KeyRotation> m_Rotations;
		std::vector<KeyScale> m_Scales;
		int m_NumPositions;
		int m_NumRotations;
		int m_NumScalings;

		DirectX::XMFLOAT4X4 m_LocalTransform;
		std::string m_Name;
		int m_ID;
	};

	class Animation
	{
	public:
		Animation(const Animation& other)
			: m_AnimationName(other.m_AnimationName),
			  m_Duration(other.m_Duration),
			  m_TicksPerSecond(other.m_TicksPerSecond),
			  m_BonesAnimationInformation(other.m_BonesAnimationInformation),
			  m_RootNode(other.m_RootNode),
			  m_BoneInfoMap(other.m_BoneInfoMap)
		{
		}

		Animation(Animation&& other) noexcept
			: m_AnimationName(std::move(other.m_AnimationName)),
			  m_Duration(other.m_Duration),
			  m_TicksPerSecond(other.m_TicksPerSecond),
			  m_BonesAnimationInformation(std::move(other.m_BonesAnimationInformation)),
			  m_RootNode(std::move(other.m_RootNode)),
			  m_BoneInfoMap(other.m_BoneInfoMap)
		{
		}

		Animation& operator=(const Animation& other)
		{
			if (this == &other)
				return *this;
			m_AnimationName = other.m_AnimationName;
			m_Duration = other.m_Duration;
			m_TicksPerSecond = other.m_TicksPerSecond;
			m_BonesAnimationInformation = other.m_BonesAnimationInformation;
			m_RootNode = other.m_RootNode;
			m_BoneInfoMap = other.m_BoneInfoMap;
			return *this;
		}

		Animation& operator=(Animation&& other) noexcept
		{
			if (this == &other)
				return *this;
			m_AnimationName = std::move(other.m_AnimationName);
			m_Duration = other.m_Duration;
			m_TicksPerSecond = other.m_TicksPerSecond;
			m_BonesAnimationInformation = std::move(other.m_BonesAnimationInformation);
			m_RootNode = std::move(other.m_RootNode);
			m_BoneInfoMap = other.m_BoneInfoMap;
			return *this;
		}


		Animation(std::string animationName, float duration, int ticksPerSecond, std::vector<BoneAnimationKeyframes> bonesAnimationInformation, std::map<std::string, BoneInfo> boneInfo, TransformNode* rootNodeOfTheModel)
		{
			m_AnimationName = animationName;
			m_RootNode = rootNodeOfTheModel;
			m_BoneInfoMap = boneInfo;
			m_BonesAnimationInformation = bonesAnimationInformation;
			m_Duration = duration;
			m_TicksPerSecond = ticksPerSecond;
		}

		~Animation()
		{
		}

		BoneAnimationKeyframes* FindBone(const std::string& name)
		{
			auto iter = std::find_if(m_BonesAnimationInformation.begin(), m_BonesAnimationInformation.end(),
			                         [&](const BoneAnimationKeyframes& Bone)
			                         {
				                         return Bone.GetBoneName() == name;
			                         }
			);
			if (iter == m_BonesAnimationInformation.end()) return nullptr;
			else return &(*iter);
		}

		inline std::string GetAnimationName() { return m_AnimationName; }
		inline float GetTicksPerSecond() { return m_TicksPerSecond; }
		inline std::vector<BoneAnimationKeyframes> GetBonesAnimationInformation() { return m_BonesAnimationInformation; }
		inline float GetDuration() { return m_Duration; }
		inline TransformNode* GetRootNode() { return m_RootNode; }
		inline void SetTicksPerSecond(int ticksPerSecond) { m_TicksPerSecond = ticksPerSecond; }
		inline const std::map<std::string, BoneInfo>& GetBoneIDMap()
		{
			return m_BoneInfoMap;
		}

	private:
		std::string m_AnimationName;
		float m_Duration;
		int m_TicksPerSecond;
		std::vector<BoneAnimationKeyframes> m_BonesAnimationInformation;
		TransformNode* m_RootNode;
		std::map<std::string, BoneInfo> m_BoneInfoMap;
	};

	class Animator
	{
	public:
		Animator()
		{
			m_CurrentTime = 0.0;
			m_CurrentAnimation = nullptr;
			m_FinalBoneMatrices.clear();
			m_FinalBoneMatrices.reserve(100);

			for (int i = 0; i < 100; i++)
			{
				DirectX::XMFLOAT4X4 identityMatrix;
				DirectX::XMStoreFloat4x4(&identityMatrix, DirectX::XMMatrixIdentity());
				m_FinalBoneMatrices.push_back(identityMatrix);
			}

		}

		void UpdateAnimation(float dt)
		{
			m_DeltaTime = dt;
			if (m_CurrentAnimation)
			{
				m_CurrentTime += m_CurrentAnimation->GetTicksPerSecond() * dt;
				m_CurrentTime = fmod(m_CurrentTime, m_CurrentAnimation->GetDuration());
				CalculateBoneTransform(m_CurrentAnimation->GetRootNode(), DirectX::XMMatrixIdentity());
			}
		}

		void PlayAnimation(Animation* pAnimation)
		{
			m_CurrentTime = 0.0;
			m_CurrentAnimation = pAnimation;
			m_FinalBoneMatrices.clear();
			m_FinalBoneMatrices.reserve(100);

			for (int i = 0; i < 100; i++)
			{
				DirectX::XMFLOAT4X4 identityMatrix;
				DirectX::XMStoreFloat4x4(&identityMatrix, DirectX::XMMatrixIdentity());
				m_FinalBoneMatrices.push_back(identityMatrix);
			}

		}

		void CalculateBoneTransform(TransformNode* node, DirectX::XMMATRIX parentTransform)
		{
			std::string nodeName = node->TransFormName();
			DirectX::XMMATRIX nodeTransform = node->GetLocalToParentSpaceTransformMatrix();

			BoneAnimationKeyframes* Bone = m_CurrentAnimation->FindBone(nodeName);

			if (Bone)
			{
				Bone->Update(m_CurrentTime);
				nodeTransform = Bone->GetLocalTransform();
			}

			DirectX::XMMATRIX globalTransformation = nodeTransform * parentTransform;

			auto boneInfoMap = m_CurrentAnimation->GetBoneIDMap();
			if (boneInfoMap.find(nodeName) != boneInfoMap.end())
			{
				int index = boneInfoMap[nodeName].BoneId;
				DirectX::XMMATRIX offset = DirectX::XMLoadFloat4x4(&boneInfoMap[nodeName].GlobalBindposeInverse);
				DirectX::XMStoreFloat4x4(&m_FinalBoneMatrices[index], offset * globalTransformation);
			}

			for (int i = 0; i < node->GetLocalChildCount(); i++)
				CalculateBoneTransform(node->GetLocalChild(i), globalTransformation);
		}

		std::vector<DirectX::XMFLOAT4X4>& GetFinalBoneMatrices()
		{
			return m_FinalBoneMatrices;
		}

	private:
		std::vector<DirectX::XMFLOAT4X4> m_FinalBoneMatrices;
		Animation* m_CurrentAnimation;
		float m_CurrentTime;
		float m_DeltaTime;

	};
}
