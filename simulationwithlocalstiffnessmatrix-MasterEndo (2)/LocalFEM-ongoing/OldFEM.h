#pragma once
//===========================================================================//
//@author Ikuto Endo
//@brief �����s���p������@�ŃV�~�����[�V�������s��Object
//===========================================================================//
#ifndef _OldFEM
#define _OldFEM

#include "Object.h"
class OldFEM : public Object {
public:
	OldFEM(std::vector<Particle*> p, ObjectData data);
	~OldFEM();

	void Update();

private:
	void Init();
	void Create_Groups();

	void Create_Group(std::vector<TetraElement*> tetra_set, int tetra_group_id);
	Eigen::Vector3f Get_Center_Grid(std::vector<Particle*> p);


	std::map<TetraElement*, bool> tetra_map;//TetraElement�����łɃO���[�v�Ɏg�p����Ă��邩�ǂ�����map
	std::vector< TetraElement* > Create_Group_Candidate();
	bool OldFEM::GroupdivideNML(TetraElement *t, int a, int b, int c);

	//debug
	unsigned int flag;

};
#endif#pragma once
#pragma once
#pragma once
