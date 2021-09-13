// naive version of collision barrier forces (ignoring heads for now), no linkedlist
double dist2;
double displX;
double displY;
double displZ;

double FXij;
double FYij;
double FZij;

const double Rnm  = 2*a;
const double Rnm2 = Rnm*Rnm;
const double chiRnm2 = 1.21*Rnm2;
double prefac;
// Collision barrier strength F^S as in paper, (21), divided by 2a.
// Must have F^S > K_B/L^2.
double FS_over_2a = 15*PI;

for (int isw=0; isw<Nsw; ++isw) {
	vec tmp = filaments[isw].getCOM();
	vec COMi(3);
	COMi << std::fmod(tmp(0),LfcmBox) << std::fmod(tmp(1),LfcmBox) << std::fmod(tmp(2),LfcmBox);


	// filament's own beads
	for(int i = 0; i<filaments[isw].length()-1; ++i) {

		// Xi = std::fmod(filaments[isw].beads[i].Xs(0),LfcmBox);
		// Yi = std::fmod(filaments[isw].beads[i].Xs(1),LfcmBox);
		// Zi = std::fmod(filaments[isw].beads[i].Xs(2),LfcmBox);

		vec XXi(3);
		XXi << std::fmod(filaments[isw].beads[i].Xs(0),LfcmBox) << std::fmod(filaments[isw].beads[i].Xs(1),LfcmBox) << std::fmod(filaments[isw].beads[i].Xs(2),LfcmBox);

		for(int j = i+2; j<filaments[isw].length(); ++j) {
			// does not include collisions with nearest neighbours

			vec XXj(3);
			XXj << std::fmod(filaments[isw].beads[j].Xs(0),LfcmBox) << std::fmod(filaments[isw].beads[j].Xs(1),LfcmBox) << std::fmod(filaments[isw].beads[j].Xs(2),LfcmBox);

			vec displ(3);
			displ = VecPeriodic(XXi,XXj);

			displX = displ(0); // Xi - filaments[isw].beads[j].Xs(0);
			displY = displ(1); // Yi - filaments[isw].beads[j].Xs(1);
			displZ = displ(2); // Zi - filaments[isw].beads[j].Xs(2);

			dist2 = (displX*displX + displY*displY + displZ*displZ);

			if(dist2 < chiRnm2) {
				prefac  =  chiRnm2 - dist2;
				prefac /= (chiRnm2 - Rnm2);
				prefac *= prefac;
				prefac *= prefac; // fourth power!
				prefac  *= FS_over_2a;

				FXij = prefac*displX;
				FYij = prefac*displY;
				FZij = prefac*displZ;

				filaments[isw].beads[i].F(0) += FXij;
				filaments[isw].beads[i].F(1) += FYij;
				filaments[isw].beads[i].F(2) += FZij;

				filaments[isw].beads[j].F(0) -= FXij;
				filaments[isw].beads[j].F(1) -= FYij;
				filaments[isw].beads[j].F(2) -= FZij;

			}
		}
	}


	// all other beads that are within range
	for (int jsw=isw+1; jsw<Nsw; ++jsw) {

		tmp = filaments[jsw].getCOM();
		vec COMj(3);
		COMj << std::fmod(tmp(0),LfcmBox) << std::fmod(tmp(1),LfcmBox) << std::fmod(tmp(2),LfcmBox);
		vec COMij = VecPeriodic(COMi,COMj);

		double comseparation = norm(COMij);

		if(comseparation < 1.5*L) {

			for(int i = 0; i<filaments[isw].length(); ++i) {

				// Xi = filaments[isw].beads[i].Xs(0);
				// Yi = filaments[isw].beads[i].Xs(1);
				// Zi = filaments[isw].beads[i].Xs(2);

				vec XXi(3);
				XXi << std::fmod(filaments[isw].beads[i].Xs(0),LfcmBox) << std::fmod(filaments[isw].beads[i].Xs(1),LfcmBox) << std::fmod(filaments[isw].beads[i].Xs(2),LfcmBox);

				for(int j = 0; j<filaments[jsw].length(); ++j) {

					// displX = Xi - filaments[jsw].beads[j].Xs(0);
					// displY = Yi - filaments[jsw].beads[j].Xs(1);
					// displZ = Zi - filaments[jsw].beads[j].Xs(2);

					vec XXj(3);
					XXj << std::fmod(filaments[isw].beads[j].Xs(0),LfcmBox) << std::fmod(filaments[isw].beads[j].Xs(1),LfcmBox) << std::fmod(filaments[isw].beads[j].Xs(2),LfcmBox);

					vec displ(3);
					displ = VecPeriodic(XXi,XXj);

					displX = displ(0); // Xi - filaments[isw].beads[j].Xs(0);
					displY = displ(1); // Yi - filaments[isw].beads[j].Xs(1);
					displZ = displ(2); // Zi - filaments[isw].beads[j].Xs(2);

					dist2 = (displX*displX + displY*displY + displZ*displZ);

					if(dist2 < chiRnm2) {
						prefac  =  chiRnm2 - dist2;
						prefac /= (chiRnm2 - Rnm2);
						prefac *= prefac;
						prefac *= prefac; // fourth power!
						prefac  *= FS_over_2a;

						FXij = prefac*displX;
						FYij = prefac*displY;
						FZij = prefac*displZ;

						filaments[isw].beads[i].F(0) += FXij;
						filaments[isw].beads[i].F(1) += FYij;
						filaments[isw].beads[i].F(2) += FZij;

						filaments[jsw].beads[j].F(0) -= FXij;
						filaments[jsw].beads[j].F(1) -= FYij;
						filaments[jsw].beads[j].F(2) -= FZij;

					}
				}
			}
		}
	}
}

return;
