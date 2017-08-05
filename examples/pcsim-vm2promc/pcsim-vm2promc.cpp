// Author: David Blyth
// pcsim-vm to ProMC converter, based on S. Chekanov's stdhep2promc converter.

#include <stdlib.h>
#include <iostream>

#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include "promc/ProMC.pb.h"
#include "promc/ProMCBook.h"
#include "promc/ProMCHeader.pb.h"

const double kEV = 100 * 1000;
const double uL = 1000;

string getEnvVar(std::string const &key) {
    char *val = getenv(key.c_str());
    return val == NULL ? std::string("") : std::string(val);
}

void readPDG(ProMCHeader *header) {
    string temp_string;
    istringstream curstring;

    string PdgTableFilename = getEnvVar("PROMC") + "/data/particle.tbl";
    if (PdgTableFilename.size() < 2) {
        cout << "**        ERROR: PROMC variable not set. Did you run source.sh"
             << "      **" << endl;
        exit(1);
    }

    ifstream fichier_a_lire(PdgTableFilename.c_str());
    if (!fichier_a_lire.good()) {
        cout << "**        ERROR: PDG Table (" << PdgTableFilename
             << ") not found! exit.                        **" << endl;
        exit(1);
        return;
    }
    // first three lines of the file are useless
    getline(fichier_a_lire, temp_string);
    getline(fichier_a_lire, temp_string);
    getline(fichier_a_lire, temp_string);
    while (getline(fichier_a_lire, temp_string)) {
        curstring.clear();  // needed when using several times istringstream::str(string)
        curstring.str(temp_string);
        long int ID;
        std::string name;
        int charge;
        float mass;
        float width;
        float lifetime;
        // ID name   chg       mass    total width   lifetime
        //  1 d      -1      0.33000     0.00000   0.00000E+00
        //  in the table, the charge is in units of e+/3
        //  the total width is in GeV
        //  the lifetime is ctau in mm
        curstring >> ID >> name >> charge >> mass >> width >> lifetime;
        ProMCHeader_ParticleData *pp = header->add_particledata();
        pp->set_id(ID);  // pi+
        pp->set_mass(mass);
        pp->set_name(name);
        pp->set_width(width);
        pp->set_lifetime(lifetime);
        pp->set_charge(charge);
        // cout << ID << " " << name << " " << mass << endl;
    }
}

void printUsage() {
    std::cerr
        << "Usage: [-d description] [-e cm_energy] [-c crosssection_in_pb] <input file> <output files>..."
        << std::endl;
}

int main(int argc, char **argv) {
    std::string description = "none";
    double cmEnergy = 0;
    double crossSection = 0;

    int opt;
    while ((opt = getopt(argc, argv, "d:e:c:h")) != -1) {
        switch (opt) {
            case 'd':
                description = optarg;
                break;
            case 'e':
                cmEnergy = atof(optarg);
                break;
            case 'c':
                crossSection = atof(optarg);
				break;
            case 'h':
                printUsage();
                exit(EXIT_SUCCESS);
            default:
                printUsage();
                exit(EXIT_FAILURE);
        }
    }
    if (argc < optind + 2) {
        printUsage();
        exit(EXIT_FAILURE);
    }

    TFile in(argv[optind]);
    TTree *tree = (TTree *)in.Get("lp_gamma_event");  // TODO: generalize
    if (!tree) return EXIT_FAILURE;

    short nParticles;
    auto array = new TClonesArray("TParticle");
    tree->SetBranchAddress("particles", &array);
    tree->SetBranchAddress("n_part", &nParticles);

    int nEntries = tree->GetEntries();
    int nFiles = argc - (optind + 1);
    char **outFiles = argv + (optind + 1);
    int entriesPerFile = nEntries / nFiles;
    std::cout << "Entries per output file: " << entriesPerFile << std::endl;

    int entry = 0;
    for (int i = 0; i < nFiles; i++) {
        ProMCBook *book = new ProMCBook(outFiles[i], "w");

        ProMCHeader header;
        header.set_momentumunit((int)kEV);
        header.set_lengthunit((int)uL);
        header.set_x1(0);
        header.set_x2(0);
        header.set_scalepdf(0);
        header.set_weight(0);
        header.set_name("pcsim-vm converted file");
        header.set_ecm(cmEnergy);
        header.set_s(cmEnergy);
        header.set_cross_section(crossSection);
        readPDG(&header);

        book->setDescription(entriesPerFile, description);
        book->setHeader(header);

        for (int j = 0; j < entriesPerFile && entry < nEntries; j++) {
            tree->GetEntry(entry);

            ProMCEvent promc;
            ProMCEvent_Event *event = promc.mutable_event();

            for (int k = 0; k < nParticles; k++) {
                auto particle = (TParticle *)(*array)[k];

                ProMCEvent_Particles *particles = promc.mutable_particles();
                particles->add_id(k);
                particles->add_pdg_id(particle->GetPDG()->PdgCode());
                particles->add_status((particle->GetWeight() == 1) ? 1 : particle->GetStatusCode());
                particles->add_px(int(particle->Px() * kEV));
                particles->add_py(int(particle->Py() * kEV));
                particles->add_pz(int(particle->Pz() * kEV));
                particles->add_energy(int(particle->Energy() * kEV));
                particles->add_mass(particle->GetMass());
                particles->add_barcode(0);
                particles->add_daughter1(particle->GetDaughter(0) + 1);
                particles->add_daughter2(particle->GetDaughter(1) + 1);
                particles->add_mother1(particle->GetMother(0) + 1);
                particles->add_mother2(particle->GetMother(1) + 1);
                particles->add_x(int(particle->Vx() * uL));
                particles->add_y(int(particle->Vy() * uL));
                particles->add_z(int(particle->Vz() * uL));
                particles->add_t(int(particle->T() * uL));
            }

            book->write(promc);
            entry++;
        }

        book->close();
    }
}
