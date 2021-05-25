#pragma once
#include <iostream>
#include <cmath>

#include "Vector.h"

namespace dem
{

	class ObjectInfo
	{
	private:

		bool is_vib_on, is_coll_on, is_visual_on;
		int id;
		double vib_amp, vib_freq;
		Vec3 post_pos, std_pos, pos;
		Vec3 h_step_velo, velo;
		Vec3 vib_unit_vec;		
		constexpr static double PI = 3.14159;

	public:

		ObjectInfo() noexcept : is_vib_on(false), is_coll_on(false), is_visual_on(false)
			, id(0), vib_amp(0.0), vib_freq(0.0)
			, post_pos(0.0, 0.0, 0.0), std_pos(0.0, 0.0, 0.0), pos(0.0, 0.0, 0.0)
			, h_step_velo(0.0, 0.0, 0.0), velo(0.0, 0.0, 0.0), vib_unit_vec(0.0, 0.0, 0.0) {}


		int GetID() const noexcept
		{
			return id; 
		}

		
		Vec3 GetPostPos() const noexcept
		{
			return post_pos;
		}


		Vec2 GetPostPosXY() const noexcept
		{
			return { post_pos.x, post_pos.y };
		}


		Vec2 GetPostPosXZ() const noexcept
		{
			return { post_pos.x, post_pos.z };
		}


		Vec2 GetPostPosYZ() const noexcept
		{
			return { post_pos.y, post_pos.z };
		}


		double GetPostPosX() const noexcept
		{
			return post_pos.x;
		}


		double GetPostPosY() const noexcept
		{
			return post_pos.y;
		}


		double GetPostPosZ() const noexcept
		{
			return post_pos.z;
		}

		
		Vec3 GetPos() const noexcept 
		{
			return pos; 
		}


		double GetPosX() const noexcept
		{
			return pos.x;
		}


		double GetPosY() const noexcept
		{
			return pos.y;
		}


		double GetPosZ() const noexcept
		{
			return pos.z;
		}


		Vec3 GetHStepVelo() const noexcept 
		{
			return h_step_velo;
		}

		
		Vec3 GetVelo() const noexcept 
		{
			return velo; 
		}	


		Vec3 GetVibUnitVec() const noexcept
		{
			return vib_unit_vec;
		}


		double GetVibAmp() const noexcept
		{
			return vib_amp;
		}

		
		double GetVibFreq() const noexcept
		{
			return vib_freq;
		}


		bool isVibOn() const noexcept
		{
			return is_vib_on; 
		}


		bool isCollOn() const noexcept
		{
			return is_coll_on; 
		}


		bool isVisualOn() const noexcept
		{
			return is_visual_on;
		}
		

		void SetID(const int v) noexcept 
		{
			id = v;
		}


		void SetVibAmp(const double v) noexcept 
		{
			vib_amp = v;
		}


		void SetVibFreq(const double v) noexcept
		{ 
			vib_freq = v;
		}


		void SetPostPos(const Vec3& v) noexcept 
		{
			post_pos = v;
		}


		void SetStdPos(const Vec3& v) noexcept
		{ 
			std_pos = v;
		}


		void SetPos(const Vec3& v) noexcept 
		{ 
			pos = v;
		}


		void SetHStepVelo(const Vec3& v) noexcept
		{
			h_step_velo = v;
		}


		void SetVelo(const Vec3& v) noexcept 
		{ 
			velo = v;
		}


		void SetVibUnitVec(const Vec3& v) noexcept
		{ 
			vib_unit_vec = v;
		}


		void SetVibFlag(const bool v) noexcept
		{
			is_vib_on = v;
		}


		void SetCollFlag(const bool v) noexcept
		{ 
			is_coll_on = v;
		}

		void SetVisualFlag(const bool v) noexcept
		{
			is_visual_on = v;
		}


		void SetHStepVelo(const double Elapsed_Time) noexcept
		{
			const double omega = 2 * PI * vib_freq;
			const Vec3 post_velo = vib_unit_vec * (vib_amp * omega * cos(omega * Elapsed_Time));
			h_step_velo = (post_velo + velo) * 0.5;
		}


		void SetPostVelo(const double Elapsed_Time) noexcept
		{
			const double omega = 2 * PI * vib_freq;
			velo = vib_unit_vec * (vib_amp * omega * cos(omega * Elapsed_Time));
		}


		void SetPostPos(const double Elapsed_Time) noexcept
		{
			pos = post_pos;
			const double omega = 2 * PI * vib_freq;
			post_pos = std_pos + (vib_unit_vec * (vib_amp * sin(omega * Elapsed_Time)));
		}

	};


	class BoxInfo : public ObjectInfo
	{
	private:

		Vec3 dim;

	public:

		BoxInfo() noexcept : ObjectInfo(), dim(0.0, 0.0, 0.0) {}
		

		Vec3 GetDim() const noexcept 
		{
			return dim;
		}


		void SetDim(const Vec3& v) noexcept 
		{
			dim = v;
		}
	};


	class HollowBoxInfo : public BoxInfo
	{
	private:

		double thickness;

	public:

		HollowBoxInfo() noexcept : BoxInfo(), thickness(0.0) {}
		

		double GetThickness() const noexcept 
		{
			return thickness; 
		}
		

		void SetThickness(const double v) noexcept 
		{
			thickness = v;
		}
		
	};


	class CylinderInfo : public ObjectInfo
	{
	private:

		double outer_radius, length;
		double rot_rad, s_plus, s_minus, c_plus, c_minus;
		Vec3 dir_unit_vec, rot_axis;

	public:

		CylinderInfo() noexcept
			: ObjectInfo(), outer_radius(0.0), length(0.0)
			, rot_rad(0.0), s_plus(0.0), s_minus(0.0), c_plus(0.0), c_minus(0.0)
			, dir_unit_vec(0.0, 0.0, 0.0), rot_axis(0.0, 0.0, 0.0) {}
		

		double GetOuterRadius() const noexcept
		{
			return outer_radius; 
		}


		double GetLength() const noexcept
		{
			return length; 
		}

		
		Vec3 GetDirUnitVec() const noexcept
		{
			return dir_unit_vec;
		}

		
		void SetOuterRadius(const double v) noexcept 
		{
			outer_radius = v;
		}


		void SetLength(const double v) noexcept 
		{
			length = v;
		}


		void SetRotRad(const double v) noexcept
		{
			rot_rad = v;
		}


		void SetSinPlus(const double v) noexcept 
		{
			s_plus = v;
		}


		void SetSinMinus(const double v) noexcept 
		{
			s_minus = v;
		}


		void SetCosPlus(const double v) noexcept
		{
			c_plus = v;
		}


		void SetCosMinus(const double v) noexcept
		{
			c_minus = v;
		}


		void SetDirUnitVec(const Vec3& v) noexcept
		{
			dir_unit_vec = v;
		}


		void SetRotAxis(const Vec3& v) noexcept 
		{
			rot_axis = v;
		}

		
		constexpr double RotatedX(const Vec3& V) const noexcept
		{
			return (V.x * c_minus) + ((rot_axis.y * V.z) - (rot_axis.z * V.y) * s_minus)
				+ ((rot_axis.x * (rot_axis.Dot(V))) * (1 - c_minus));
		}
		

		constexpr Vec3 Rotated(const Vec3& V) const noexcept
		{
			return{ (V * c_minus) + (rot_axis.Cross(V) * s_minus)
				+ ((rot_axis * (rot_axis.Dot(V))) * (1 - c_minus)) };
		}
		

		constexpr Vec3 RotatedReturn(const Vec3& V) const noexcept
		{
			return{ (V * c_plus) + (rot_axis.Cross(V) * s_plus)
				+ ((rot_axis * (rot_axis.Dot(V))) * (1 - c_plus)) };
		}

	};


	class HollowCylinderInfo : public CylinderInfo
	{
	private:

		double inner_radius;

	public:

		HollowCylinderInfo() noexcept : CylinderInfo(), inner_radius(0.0) {}
		

		double GetInnerRadius() const noexcept 
		{
			return inner_radius; 
		}


		void SetInnerRadius(const double v) noexcept
		{
			inner_radius = v; 
		}

	};


	class SphereInfo : public ObjectInfo
	{
	private:

		double outer_radius;

	public:

		SphereInfo() noexcept : ObjectInfo(), outer_radius(0.0) {}
		

		double GetOuterRadius() const noexcept 
		{
			return outer_radius; 
		}


		void SetOuterRadius(const double v) noexcept
		{
			outer_radius = v;
		}
	};


	class HollowSphereInfo : public SphereInfo
	{
	private:

		double inner_radius;

	public:

		HollowSphereInfo() noexcept	: SphereInfo(), inner_radius(0.0) {}

		
		double GetInnerRadius() const noexcept
		{
			return inner_radius; 
		}
		

		void SetInnerRadius(const double v) noexcept
		{
			inner_radius = v; 
		}

	};

}