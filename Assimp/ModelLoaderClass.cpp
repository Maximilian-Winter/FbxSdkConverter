#include "ModelLoaderClass.h"



SurfaceMaterialData::SurfaceMaterialData()
{
	matName = "default";
	AlbedoMapFilename = "defaultAlbedo.dds";
	RoughnessMapFilename = "defaultRoughness.dds";
	MetalnessMapFilename = "defaultMetalness.dds";
	f0 = 0.06;
	RenderForward = false;
	HasNormalMap = false;
	useAlpaChannelTransparency = false;
	is_transparent = false;
	materialIndex = 0;
}

Vertex::Vertex()
{
	Position.x = 0.0;
	Position.y = 0.0;
	Position.z = 0.0;
	Normal.x = 0.0;
	Normal.y = 0.0;
	Normal.z = 0.0;
	TexCoords.x = 0.0;
	TexCoords.y = 0.0;
}

Vertex::Vertex(DirectX::XMFLOAT3 VertexPosition, DirectX::XMFLOAT3 VertexNormal, DirectX::XMFLOAT2 VertexTexCoords)
{
	Position = VertexPosition;
	Normal = VertexNormal;
	TexCoords = VertexTexCoords;
}

AnimatedVertex::AnimatedVertex()
{
	Position.x = 0.0;
	Position.y = 0.0;
	Position.z = 0.0;
	Normal.x = 0.0;
	Normal.y = 0.0;
	Normal.z = 0.0;
	TexCoords.x = 0.0;
	TexCoords.y = 0.0;

}

AnimatedVertex::AnimatedVertex(DirectX::XMVECTOR VertexPosition, DirectX::XMVECTOR VertexNormal,
	DirectX::XMVECTOR VertexTexCoords, std::list<BoneIdBoneWeight> boneIDBoneWeight)
{
	DirectX::XMStoreFloat3(&Position, VertexPosition);
	DirectX::XMStoreFloat3(&Normal, VertexNormal);
	DirectX::XMStoreFloat2(&TexCoords, VertexTexCoords);
	BoneIDBoneWeight = boneIDBoneWeight;
}

TriangleVertexIndicies::TriangleVertexIndicies()
{
	Indices[0] = 0;
	Indices[1] = 0;
	Indices[2] = 0;
}

TriangleVertexIndicies::TriangleVertexIndicies(int indices[3])
{
	Indices[0] = indices[0];
	Indices[1] = indices[1];
	Indices[2] = indices[2];
}

Mesh::Mesh(std::vector<Vertex> vertices, std::vector<TriangleVertexIndicies> indices, SubsetSurfaceData subsetData)
{

	this->Vertices = vertices;
	this->Indices = indices;
	this->SubsetData = subsetData;
}

Mesh::Mesh(std::string meshName, std::vector<Vertex> vertices, std::vector<TriangleVertexIndicies> indices,
	SubsetSurfaceData subsetData)
{
	this->MeshName = meshName;
	this->Vertices = vertices;
	this->Indices = indices;
	this->SubsetData = subsetData;
}

Vertex Mesh::GetVertexData(int VertexIndex)
{
	if(VertexIndex > Vertices.size())
	{
		std::cout << "TEST" << std::endl;
	}
	return Vertices[VertexIndex];
}

int Mesh::GetNumberOfVertices()
{
	return Vertices.size();
}

TriangleVertexIndicies Mesh::GetTriangleIndexData(int TriangleIndex)
{
	return Indices[TriangleIndex];
}

int Mesh::GetNumberOfTriangleIndicies()
{
	return Indices.size();
}


AnimatedMesh::AnimatedMesh(std::vector<AnimatedVertex> vertices, std::vector<TriangleVertexIndicies> indices, std::vector<SubsetSurfaceData> subsetData)
{

	this->Vertices = vertices;
	this->Indices = indices;
	this->SubsetData = subsetData;
}

AnimatedMesh::AnimatedMesh(std::string meshName, std::vector<AnimatedVertex> vertices,
	std::vector<TriangleVertexIndicies> indices, std::vector<SubsetSurfaceData> subsetData)
{
	this->MeshName = meshName;
	this->Vertices = vertices;
	this->Indices = indices;
	this->SubsetData = subsetData;
}

AnimatedVertex AnimatedMesh::GetVertexData(int VertexIndex)
{
	if (VertexIndex > Vertices.size())
	{
		std::cout << "TEST" << std::endl;
	}
	return Vertices[VertexIndex];
}

int AnimatedMesh::GetNumberOfVertices()
{
	return Vertices.size();
}

TriangleVertexIndicies AnimatedMesh::GetTriangleIndexData(int TriangleIndex)
{
	return Indices[TriangleIndex];
}

int AnimatedMesh::GetNumberOfTriangleIndicies()
{
	return Indices.size();
}

void StringUtilities::RemoveCharsFromString(std::string& str, char charToRemove)
{
	str.erase(std::ranges::remove(str, charToRemove).begin(), str.end());
}

