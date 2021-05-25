#include "Object.h"

namespace dem
{

	void Object::ShowValue() noexcept
	{
		ShowValueFloorBase();
		ShowValueFieldBase();
		ShowValueBoundaryBase();
		ShowValueBox();
		ShowValueHollowBox();
		ShowValueCylinder();
		ShowValueHollowCylinder();
		ShowValueSphere();
		ShowValueHollowSphere();
	}


	void Object::ShowValueFloorBase() noexcept
	{
		for (const auto& o : floor_base)
		{
			std::cout << std::endl;
			std::cout << "floor_base : " << std::endl;
			ShowValue(o);
		}
	}


	void Object::ShowValueFieldBase() noexcept
	{

		for (const auto& o : field_base)
		{
			std::cout << std::endl;
			std::cout << "field_base : " << std::endl;
			ShowValue(o);
		}
	}


	void Object::ShowValueBoundaryBase() noexcept
	{

		for (const auto& o : boundary_base)
		{
			std::cout << std::endl;
			std::cout << "boundary_base : " << std::endl;
			ShowValue(o);
		}
	}


	void Object::ShowValueBox() noexcept
	{
		for (const auto& o : box)
		{
			std::cout << std::endl;
			std::cout << "box : " << std::endl;
			ShowValue(o);
			std::cout << "Dim (" << o.GetDim() << ")" << std::endl;
		}
	}


	void Object::ShowValueHollowBox() noexcept
	{
		for (const auto& o : hollow_box)
		{
			std::cout << std::endl;
			std::cout << "hollow_box : " << std::endl;
			ShowValue(o);
			std::cout << "Dim (" << o.GetDim() << ")" << std::endl;
			std::cout << "Thickness (" << o.GetThickness() << ")" << std::endl;
		}
	}


	void Object::ShowValueCylinder() noexcept
	{
		for (const auto& o : cylinder)
		{
			std::cout << std::endl;
			std::cout << "cylinder : " << std::endl;
			ShowValue(o);
			std::cout << "Outer Radius (" << o.GetOuterRadius() << ")" << std::endl;
			std::cout << "Length (" << o.GetLength() << ")" << std::endl;
			std::cout << "Dir_Unit_Vec (" << o.GetDirUnitVec() << ")" << std::endl;
		}
	}


	void Object::ShowValueHollowCylinder() noexcept
	{
		for (const auto& o : hollow_cylinder)
		{
			std::cout << std::endl;
			std::cout << "hollow_cylinder : " << std::endl;
			ShowValue(o);
			std::cout << "Inner Radius (" << o.GetInnerRadius() << ")" << std::endl;
			std::cout << "Outer Radius (" << o.GetOuterRadius() << ")" << std::endl;
			std::cout << "Length (" << o.GetLength() << ")" << std::endl;
			std::cout << "Dir_Unit_Vec (" << o.GetDirUnitVec() << ")" << std::endl;
		}
	}


	void Object::ShowValueSphere() noexcept
	{
		for (const auto& o : sphere)
		{
			std::cout << std::endl;
			std::cout << "sphere : " << std::endl;
			ShowValue(o);
			std::cout << "Outer Radius (" << o.GetOuterRadius() << ")" << std::endl;			
		}
	}


	void Object::ShowValueHollowSphere() noexcept
	{
		for (const auto& o : hollow_sphere)
		{
			std::cout << std::endl;
			std::cout << "hollow_sphere : " << std::endl;
			ShowValue(o);
			std::cout << "Inner Radius (" << o.GetInnerRadius() << ")" << std::endl;
			std::cout << "Outer Radius (" << o.GetOuterRadius() << ")" << std::endl;
		}
	}


	void Object::SetHStepVelo(const double Elapsed_Time) noexcept
	{
		for (auto& o : floor_base)	   if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : field_base)     if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : boundary_base)  if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : box)			   if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : hollow_box)     if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : cylinder)       if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : hollow_cylinder)if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : sphere)         if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
		for (auto& o : hollow_sphere)  if (o.isVibOn()) o.SetHStepVelo(Elapsed_Time);
	}


	void Object::SetPostVelo(const double Elapsed_Time) noexcept
	{
		for (auto& o : floor_base)	   if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : field_base)     if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : boundary_base)  if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : box)			   if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : hollow_box)     if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : cylinder)       if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : hollow_cylinder)if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : sphere)         if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
		for (auto& o : hollow_sphere)  if (o.isVibOn()) o.SetPostVelo(Elapsed_Time);
	}


	void Object::SetPostPos(const double Elapsed_Time) noexcept
	{
		for (auto& o : floor_base)     if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : field_base)     if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : boundary_base)  if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : box)            if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : hollow_box)     if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : cylinder)       if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : hollow_cylinder)if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : sphere)         if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
		for (auto& o : hollow_sphere)  if (o.isVibOn()) o.SetPostPos(Elapsed_Time);
	}
		

	void Object::LoadObj(const std::vector<std::string>& Obj_Filename) noexcept
	{
		//open the object file
		std::cout << "loading object file..." << std::endl;

		for (const auto& o : Obj_Filename)
		{
			std::cout << "open [" << o << "]" << std::endl;
			std::ifstream fin(o);
			if (!fin)	data_io.ErrorInput(o);
			
			std::string line, obj_type;
			std::vector<std::string> obj_info;
			bool obj_on = 0, first_line = 1;

			while (!fin.eof() && getline(fin, line))
			{
				std::vector<std::string> key;
				data_io.SplitLine(line, key);
				if (key.empty()) continue;

				if ((key[0] == "Object_End") && (obj_on ))
				{
					obj_on = 0;
					CreateObj(obj_type, obj_info);
				}

				if (obj_on)
				{
					if (first_line) 
					{
						obj_type = key[0];
						first_line = 0; 
					}
					else { obj_info.push_back(line); }
				}

				if (key[0] == "Object")
				{
					obj_on = 1; first_line = 1;
					obj_info.clear();
				}
			}
			fin.close();
		}

		// boundary_base
		if (!GetBoundaryBaseSize())
		{
			ObjectInfo object;
			boundary_base.push_back(object);
		}

		// boundary_base
		if (!GetFieldBaseSize())
		{
			ObjectInfo object;
			field_base.push_back(object);
		}

		ShowValue();
		std::cout << "finish loading object data" << std::endl;
	}


	void Object::CreateObj(const std::string& Obj_Type,
		const std::vector<std::string>& Obj_Info) noexcept
	{
		std::map<std::string, std::string> parameter;
		SetObjParameter(parameter, Obj_Info);
		const int parameter_size = static_cast<int>(parameter.size());

		if (Obj_Type == "Floor")
		{
			if ((parameter_size == 11) && (GetFloorBaseSize() == 0))
			{
				ObjectInfo object;
				object.SetID(0);
				SetBasicParameter(object, parameter);	
				SetCollFlagVisualFlag(object, parameter);
				floor_base.push_back(object);
			}
			else { data_io.ErrorInput("Floor"); }
		}

		if (Obj_Type == "Field_Base")
		{
			if ((parameter_size == 9) && (GetFieldBaseSize() == 0))
			{
				ObjectInfo object;
				object.SetID(0);
				SetBasicParameter(object, parameter);
				field_base.push_back(object);
			}
			else { data_io.ErrorInput("Field_Base"); }
		}

		if (Obj_Type == "Boundary_Base")
		{
			if ((parameter_size == 9) && (GetBoundaryBaseSize() == 0))
			{
				ObjectInfo object;
				object.SetID(0);
				SetBasicParameter(object, parameter);
				boundary_base.push_back(object);
			}
			else { data_io.ErrorInput("Boundary_Base"); }
		}

		if (Obj_Type == "Box")
		{
			if (parameter_size == 14)
			{
				BoxInfo object;
				object.SetID(GetBoxSize());
				SetBasicParameter(object, parameter);
				SetCollFlagVisualFlag(object, parameter);
				SetBoxParameter(object, parameter);
				box.push_back(object);
			}
			else { data_io.ErrorInput("Box"); }
		}

		if (Obj_Type == "Hollow_Box")
		{
			if (parameter_size == 15)
			{
				HollowBoxInfo object;
				object.SetID(GetHollowBoxSize());
				SetBasicParameter(object, parameter);
				SetCollFlagVisualFlag(object, parameter);
				SetBoxParameter(object, parameter);
				object.SetThickness(data_io.StringToDouble(parameter["thickness"]));
				hollow_box.push_back(object);
			}
			else { data_io.ErrorInput("Hollow_Box"); }
		}

		if (Obj_Type == "Cylinder")
		{
			if (parameter_size == 16)
			{
				CylinderInfo object;
				object.SetID(GetCylinderSize());
				SetBasicParameter(object, parameter);
				SetCollFlagVisualFlag(object, parameter);
				SetCylinderParameter(object, parameter);
				SetOuterRadius(object, parameter);
				cylinder.push_back(object);
			}
			else { data_io.ErrorInput("Cylinder"); }
		}

		if (Obj_Type == "Hollow_Cylinder")
		{
			if (parameter_size == 17)
			{
				HollowCylinderInfo object;
				object.SetID(GetHollowCylinderSize());
				SetBasicParameter(object, parameter);
				SetCollFlagVisualFlag(object, parameter);
				SetCylinderParameter(object, parameter);
				SetInnerRadius(object, parameter);
				SetOuterRadius(object, parameter);
				hollow_cylinder.push_back(object);
			}
			else { data_io.ErrorInput("Hollow_Cylinder"); }
		}

		if (Obj_Type == "Sphere")
		{
			if (parameter_size == 12)
			{
				SphereInfo object;
				object.SetID(GetSphereSize());
				SetBasicParameter(object, parameter);
				SetCollFlagVisualFlag(object, parameter);
				SetOuterRadius(object, parameter);
				sphere.push_back(object);
			}
			else { data_io.ErrorInput("Sphere"); }
		}

		if (Obj_Type == "Hollow_Sphere")
		{
			if (parameter_size == 13)
			{
				HollowSphereInfo object;
				object.SetID(GetHollowSphereSize());
				SetBasicParameter(object, parameter);
				SetCollFlagVisualFlag(object, parameter);
				SetInnerRadius(object, parameter);
				SetOuterRadius(object, parameter);
				hollow_sphere.push_back(object);
			}
			else { data_io.ErrorInput("Hollow_Sphere"); }
		}
	}


	void Object::SetObjParameter(std::map<std::string, std::string>& Parameter, 
		const std::vector<std::string>& Obj_Info) noexcept
	{
		for (const auto& i : Obj_Info)
		{
			std::vector<std::string> key;
			data_io.SplitLine(i, key);
			const int key_size = static_cast<int>(key.size());

			if (key[0] == "Position")
			{
				if (key_size == 4)
				{
					Parameter["pos_x"] = key[1];
					Parameter["pos_y"] = key[2];
					Parameter["pos_z"] = key[3];
				}
				else { data_io.ErrorInput("Position"); }
			}

			if (key[0] == "Outer_Radius")
			{
				if (key_size == 2) { Parameter["outer_radius"] = key[1]; }
				else { data_io.ErrorInput("Outer_Radius"); }
			}

			if (key[0] == "Inner_Radius")
			{
				if (key_size == 2) { Parameter["inner_radius"] = key[1]; }
				else { data_io.ErrorInput("Inner_Radius"); }
			}

			if (key[0] == "Length")
			{
				if (key_size == 2) { Parameter["length"] = key[1]; }
				else { data_io.ErrorInput("Length"); }
			}

			if (key[0] == "Thickness")
			{
				if (key_size == 2) { Parameter["thickness"] = key[1]; }
				else { data_io.ErrorInput("Thickness"); }
			}

			if (key[0] == "Dimention")
			{
				if (key_size == 4)
				{
					Parameter["dim_x"] = key[1];
					Parameter["dim_y"] = key[2];
					Parameter["dim_z"] = key[3];
				}
				else { data_io.ErrorInput("Dimention"); }
			}

			if (key[0] == "Direction_Unit_Vector")
			{
				if (key_size == 4)
				{
					Parameter["dir_x"] = key[1];
					Parameter["dir_y"] = key[2];
					Parameter["dir_z"] = key[3];
				}
				else { data_io.ErrorInput("Direction_Unit_Vector"); }
			}

			if (key[0] == "Collision_On")
			{
				if (key_size == 2) { Parameter["collision_on"] = key[1]; }
				else { data_io.ErrorInput("Collision_On"); }
			}
			
			if (key[0] == "Vibration_On")
			{
				if (key_size == 2) { Parameter["vibration_on"] = key[1]; }
				else { data_io.ErrorInput("Vibration_On"); }
			}

			if (key[0] == "Visualization_On")
			{
				if (key_size == 2) { Parameter["visualization_on"] = key[1]; }
				else { data_io.ErrorInput("Visualization_On"); }
			}

			if (key[0] == "Vibration_Frequency")
			{
				if (key_size == 2) { Parameter["frequency"] = key[1]; }
				else { data_io.ErrorInput("Vibration_Frequency"); }
			}

			if (key[0] == "Vibration_Amplitude")
			{
				if (key_size == 2) { Parameter["amplitude"] = key[1]; }
				else { data_io.ErrorInput("Vibration_Amplitude"); }
			}

			if (key[0] == "Vibration_Unit_Vector")
			{
				if (key_size == 4)
				{
					Parameter["vib_x"] = key[1];
					Parameter["vib_y"] = key[2];
					Parameter["vib_z"] = key[3];
				}
				else { data_io.ErrorInput("Vibration_Unit_Vector"); }
			}
		}
	}

}

