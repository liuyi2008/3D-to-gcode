#include <algorithm>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#ifndef DOXY_IGNORE_THIS

struct MyTraits : public OpenMesh::DefaultTraits
{
	HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
	VertexTraits
	{
		int some_additional_index;
	};
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

bool operator <(const MyMesh::VertexIter & a1, const MyMesh::VertexIter & a2)
{
	//MyMesh::Point p1 = mesh.point(a1);
	//MyMesh::Point p2 = mesh.point(a2);
	//float p1z = p1[3];

	//return p1[3] > p2[3];
	return true;
}

template <class Mesh> class SmootherT
{
public:
	typedef typename Mesh::Point            cog_t;
	typedef OpenMesh::VPropHandleT< cog_t > Property_cog;
public:
	// construct with a given mesh
	explicit SmootherT(Mesh& _mesh)
		: mesh_(_mesh)
	{
		mesh_.add_property(cog_);
	}
	~SmootherT()
	{
		mesh_.remove_property(cog_);
	}
	// smooth mesh _iterations times
	void smooth(unsigned int _iterations)
	{
		for (unsigned int i = 0; i < _iterations; ++i)
		{
			std::for_each(mesh_.vertices_begin(), mesh_.vertices_end(), ComputeCOG(mesh_, cog_));
			std::for_each(mesh_.vertices_begin(),
				mesh_.vertices_end(),
				SetCOG(mesh_, cog_));

			sort(mesh.vertices_begin(), mesh.vertices_end());
		}
	}

	void sortVertZ()
	{
		//sort(mesh.vertices_begin(), mesh.vertices_end());
		//sort(mesh.vertices_begin(), mesh.vertices_end(), GreaterVertZ);
	}
private:
	//--- private classes ---
	class ComputeCOG
	{
	public:
		ComputeCOG(Mesh& _mesh, Property_cog& _cog)
			: mesh_(_mesh), cog_(_cog)
		{}
		void operator()(const typename Mesh::VertexHandle& _vh)
		{
			typename Mesh::VertexVertexIter  vv_it;
			typename Mesh::Scalar            valence(0.0);

			mesh_.property(cog_, _vh) = typename Mesh::Point(0.0, 0.0, 0.0);
			for (vv_it = mesh_.vv_iter(_vh); vv_it.is_valid(); ++vv_it)
			{
				mesh_.property(cog_, _vh) += mesh_.point(*vv_it);
				++valence;
			}
			mesh_.property(cog_, _vh) /= valence;
		}
	private:
		Mesh&         mesh_;
		Property_cog& cog_;
	};
	class SetCOG
	{
	public:
		SetCOG(Mesh& _mesh, Property_cog& _cog)
			: mesh_(_mesh), cog_(_cog)
		{}
		void operator()(const typename Mesh::VertexHandle& _vh)
		{
			if (!mesh_.is_boundary(_vh))
				mesh_.set_point(_vh, mesh_.property(cog_, _vh));
		}
	private:
		Mesh&         mesh_;
		Property_cog& cog_;
	};

	//--- private elements ---
	Mesh&        mesh_;
	Property_cog cog_;
};
#endif
