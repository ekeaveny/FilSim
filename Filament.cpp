#include "multi_filament_header.hpp"

int Filament::length() const {
	return myNworm;
}

// set up in random initial orientation given head position
Filament& Filament::initialFilamentSetupIsoRand(const vec& HeadPos,
                                                const int forceArrayLengthFilled,
                                                const int stateArrayLengthFilled,
                                                double phase_in) {
	vec X(3);
	vec tangent(3);

	Phase = phase_in;

	beads[0].initialBeadRandom(HeadPos,ah);
	beads[0].qs.assignTangentVec(tangent);

	for(int n = 1; n<myNworm; ++n) {
		X = HeadPos + (n*dL)*tangent;
		beads[n].initialBeadSetupQuat(X,beads[0].qs,a);
	}
	this->forceLengthFilled = forceArrayLengthFilled;
	this->stateLengthFilled = stateArrayLengthFilled;

	return *this;
}

// set up all beads in an initially straight filament given head's position
// tangent vector and normal vector.
Filament& Filament::initialFilamentSetup(const vec& HeadPos,
                                         const vec& tangent, const vec& normal,
                                         const int forceArrayLengthFilled,
                                         const int stateArrayLengthFilled) {
	vec X(3);
	for(int n = 0; n<myNworm; ++n) {
		X = HeadPos + (n*dL)*tangent; //position of n-th bead of swimmer
		beads[n].initialBeadSetup(X,tangent,normal,a);
	}
	this->forceLengthFilled = forceArrayLengthFilled;
	this->stateLengthFilled = stateArrayLengthFilled;
	return *this;
}

// random initial orientation with the centre at the spawn point
Filament& Filament::initialFilamentSetupIsoRandCentred(const vec& HeadPos,
                                                       const int forceArrayLengthFilled,
                                                       const int stateArrayLengthFilled) {
	vec X(3);
	vec tangent(3);

	beads[0].initialBeadRandom(HeadPos,ah);
	beads[0].qs.assignTangentVec(tangent);

	for(int n = 0; n<myNworm; ++n) {
		//position of n-th bead of swimmer, subtracting off half the length
		X = HeadPos + (n*dL)*tangent - myNworm*dL/2*tangent;
		beads[n].initialBeadSetupQuat(X,beads[0].qs,a);
	}

	this->forceLengthFilled = forceArrayLengthFilled;
	this->stateLengthFilled = stateArrayLengthFilled;

	return *this;
}

// random initial orientation with the centre at the spawn point
Filament& Filament::initialFilamentSetupBetweenTwoPoints(const vec& node1,
                                                         const vec& node2,
                                                         const int forceArrayLengthFilled,
                                                         const int stateArrayLengthFilled){

	vec X(3);
	vec tangent(3), normal(3);
	int num_beads_required;
	vec CoM(3), unit_vector(3);

	tie(num_beads_required, CoM, unit_vector) =
		beads_between_nodes(node1, node2);

	// Number of beads should already be set by the filament constructor

	vec binormal = {1,0,0}; // Arbitrary, I think
	tangent = unit_vector;
	normal = cross(tangent,binormal);
	normal = normal/norm(normal);

	for(int n = 0; n<myNworm; ++n) {
		//position of n-th bead of swimmer, subtracting off half the length
		X = CoM + (n*dL)*tangent - myNworm*dL/2*tangent;
		beads[n].initialBeadSetup(X,tangent,normal,a);
	}

	this->forceLengthFilled = forceArrayLengthFilled;
	this->stateLengthFilled = stateArrayLengthFilled;

	return *this;
}

// Check overlap with all other filaments in array otherFilaments[]
// until NumberOthers array element
bool Filament::checkOverlapWith(const std::vector<Filament>& otherFilaments,
                                const int NumberOthers,
                                float sep_factor) const {
	// loop over all other filaments' beads
	for(int n = 0; n<NumberOthers; ++n) {
		//check centre of mass distance
		if(norm(this->getCOM() - otherFilaments[n].getCOM()) < 2.5*L) {
			for(int i = 0; i<otherFilaments[n].myNworm; ++i) {
				for(int j = 0; j<myNworm; ++j) {
					if(norm(beads[j].X - otherFilaments[n].beads[i].X) <
					   sep_factor*(beads[j].radius + otherFilaments[n].beads[i].radius)) { // this is a big initial separation!
						return true;
					}
				}
			}
		}
	}
	return false;
}

// Check overlap with all other filaments in array otherFilaments[]
// until NumberOthers array element. Takes into account periodic boundary
// conditions. Not very efficiently written, as called only once.
bool Filament::checkOverlapWithPeriodic(const std::vector<Filament>& otherFilaments,
                                        const int NumberOthers,
                                        float sep_factor)  {
	// sep_factor: default value = 2

	vec com1(3);
	vec com2(3);

	for(int n = 0; n<NumberOthers; ++n) {
		//check centre of mass distance
		com1 = this->getCOM();
		com2 = otherFilaments[n].getCOM();

		if(norm(VecPeriodicX(com1, com2)) < 2.5*L) {
			for(int i = 0; i<otherFilaments[n].myNworm; ++i) {
				for(int j = 0; j<myNworm; ++j) {
					double distance_apart = norm(VecPeriodicX(beads[j].X,otherFilaments[n].beads[i].X));
					double separation_allowed = sep_factor*(beads[j].radius + otherFilaments[n].beads[i].radius);
					if (distance_apart < separation_allowed) {
						cout << endl << " Fil " << n << " bead " << i << " overlaps with proposed bead " << j << endl;
						return true;
					}
				}
			}
		}
	}
	return false;
}
