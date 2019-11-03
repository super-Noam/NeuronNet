#include <gtest/gtest.h>
#include <algorithm>
#include "random.h"
#include "simulation.h"

RandomNumbers *_RNG = new RandomNumbers(101301091);
Network net;
Neuron n1, n2;
size_t nlinks = 1000;
double noise = 10.;

TEST(simulationTest, parse) {
    Simulation s1(10, 1);
    s1.parse_types("");
    EXPECT_EQ(2, s1.size_type("FS"));
    EXPECT_EQ(8, s1.size_type("RS"));
    EXPECT_EQ(0, s1.size_type("XX"));
    s1.parse_types("IB:0.5, CH:0.2,TC:0.1");
    EXPECT_EQ(2, s1.size_type("RS"));    
    EXPECT_EQ(2, s1.size_type("CH"));
    EXPECT_EQ(5, s1.size_type("IB"));
}

TEST(neuronTest, initialize) {
    n1.set_default_params("RS");
    n2.set_default_params("FS");
    EXPECT_FALSE(n1.is_inhibitory());
    EXPECT_TRUE(n2.is_inhibitory());
    EXPECT_TRUE(n2.is_type("FS"));
    EXPECT_FALSE(n2.is_type("RS"));
    EXPECT_EQ("RS\t0.02\t0.2\t-65\t8\t0", n1.formatted_params());
    EXPECT_EQ("FS\t0.1\t0.2\t-65\t2\t1",  n2.formatted_params());
}

TEST(neuronTest, step) {
    n1.input(10.0);
    n1.step();
    EXPECT_NEAR(-58.105, n1.potential(), 1e-4);
    EXPECT_NEAR(-12.9724, n1.recovery(), 1e-4);
    EXPECT_FALSE(n1.firing());
    for (int n=0; n<3; n++) {
        n1.input(10.0);
        n1.step();
    }
    EXPECT_TRUE(n1.firing());
}

TEST(networkTest, initialize) {
    net.resize(nlinks);
    EXPECT_EQ(nlinks, net.size());
    double mean = 0, sdv = 0;
    std::vector<double> allv(net.recoveries());
    for (auto I : allv) {
        mean += 0.001*I;
        sdv  += 0.001*I*I;
    }
    sdv -= mean*mean;	
    sdv = std::sqrt(3*sdv);
// --- recovery potential is -65 (rest value) X 0.2 (default b) = -13 for excitatory
// --- and -13 * (1+_BVAR_/2) for inhibitory
    EXPECT_NEAR(0.2*_REST_VAL_*(1+0.1*_BVAR_), mean, 1e-1);
    EXPECT_NEAR(-0.2*_REST_VAL_*_BVAR_*sqrt(0.194), sdv, 1.5e-1);
}

TEST(networkTest, connect) {
    bool trylink = net.add_link(0, 0, 10);
    EXPECT_FALSE(trylink);
    
    trylink = net.add_link(0, nlinks+1, .5);
    EXPECT_FALSE(trylink);
    
    size_t excit_idx = 0;
    for (; excit_idx<net.size(); excit_idx++)	
        if (!net.neuron(excit_idx).is_inhibitory())
            break;

    size_t nlink = 3;
    double stren = 6.;
    std::pair<size_t, double> expec{nlink, -2*stren*nlink};
    
// --- add 3 links from inhibitors to same excitatory

    for (size_t inhib_idx=0; inhib_idx<net.size() && nlink>0; inhib_idx++) 
        if ((net.neuron(inhib_idx).is_inhibitory())
            && net.add_link(excit_idx, inhib_idx, stren))
            nlink--;

    EXPECT_EQ(expec, net.degree(excit_idx));

    std::vector<double> noisev(nlinks, noise);
    size_t inhib1 = net.neighbors(excit_idx).front().first;		
    int ifirs(0), efirs(0);
// --- should fire once within 10 time steps
    for (size_t t=0; t<10; t++) {
        net.step(noisev);
        ifirs += (int)net.neuron(inhib1).firing();
        efirs += (int)net.neuron(excit_idx).firing();
    }
    EXPECT_EQ(1, ifirs);
    
    EXPECT_EQ(1, efirs);
    
// --- input to excitatory should be noise
// --- input to inhibitory should be .4*noise
    EXPECT_DOUBLE_EQ(noise,    net.neuron(excit_idx).input());
    EXPECT_DOUBLE_EQ(.4*noise, net.neuron(inhib1).input());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
