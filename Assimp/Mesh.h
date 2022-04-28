#pragma once

#include <map>

//Includes
#include <DirectXMath.h>
#include <string>
#include <vector>
#include <iostream>

#include <btBulletDynamicsCommon.h>
#include <BulletCollision\CollisionShapes\btShapeHull.h>
#include <BulletCollision\CollisionShapes\btConvexHullShape.h>
#include <Extras\ConvexDecomposition\ConvexDecomposition.h>
#include <Extras\HACD\hacdHACD.h>

#include <LinearMath\btSerializer.h>
#include <Extras\Serialize\BulletWorldImporter\btBulletWorldImporter.h>

//#include <VHACD.h>

struct SurfaceMaterialData
{
	SurfaceMaterialData();
	std::string matName;
	std::string AlbedoMapFilename;
	std::string RoughnessMapFilename;
	std::string MetalnessMapFilename;
	std::string NormalMapFilename;
	int materialIndex;
	float f0;
	bool RenderForward;
	bool useAlpaChannelTransparency;
	bool is_transparent;
	bool HasNormalMap;
};

struct Subset
{
	Subset() :
		materialName("Default"),
		start(0),
		drawAmount(0)
	{
	}

	Subset(std::string materialName, int start, int drawAmount)
		:materialName(materialName), start(start), drawAmount(drawAmount)
	{
	}
	std::string materialName;
	int start;
	int drawAmount;
};


struct Vertex
{
	Vertex();

	Vertex(DirectX::XMFLOAT3 VertexPosition, DirectX::XMFLOAT3 VertexNormal, DirectX::XMFLOAT2 VertexTexCoords);

	DirectX::XMFLOAT3 Position{};
	DirectX::XMFLOAT3 Normal{};
	DirectX::XMFLOAT2 TexCoords{};
};


struct SubsetSurfaceData
{
	SubsetSurfaceData() = default;

	SubsetSurfaceData(const SubsetSurfaceData& other)
		: SurfaceMaterial(other.SurfaceMaterial),
		  Start(other.Start),
		  DrawAmount(other.DrawAmount)
	{
	}

	SubsetSurfaceData(SubsetSurfaceData&& other) noexcept
		: SurfaceMaterial(std::move(other.SurfaceMaterial)),
		  Start(other.Start),
		  DrawAmount(other.DrawAmount)
	{
	}

	SubsetSurfaceData& operator=(const SubsetSurfaceData& other)
	{
		if (this == &other)
			return *this;
		SurfaceMaterial = other.SurfaceMaterial;
		Start = other.Start;
		DrawAmount = other.DrawAmount;
		return *this;
	}

	SubsetSurfaceData& operator=(SubsetSurfaceData&& other) noexcept
	{
		if (this == &other)
			return *this;
		SurfaceMaterial = std::move(other.SurfaceMaterial);
		Start = other.Start;
		DrawAmount = other.DrawAmount;
		return *this;
	}

	SurfaceMaterialData SurfaceMaterial;
	int Start;
	int DrawAmount;
};

struct TriangleVertexIndicies
{
	TriangleVertexIndicies();

	TriangleVertexIndicies(int indices[3]);

	int Indices[3];
};

struct Mesh
{
	Mesh(std::vector<Vertex> vertices, std::vector<TriangleVertexIndicies> indices, SubsetSurfaceData subsetData);

	Mesh(std::string meshName, std::vector<Vertex> vertices, std::vector<TriangleVertexIndicies> indices, SubsetSurfaceData subsetData);

	Vertex GetVertexData(int VertexIndex);
	int GetNumberOfVertices();

	TriangleVertexIndicies GetTriangleIndexData(int TriangleIndex);
	int GetNumberOfTriangleIndicies();

	std::string MeshName;

	std::vector<Vertex> Vertices;
	std::vector<TriangleVertexIndicies> Indices;
	SubsetSurfaceData SubsetData;
};


