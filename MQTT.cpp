#pragma once
#include "MQTT.hpp"
#include "MQTTScheduler.hpp"
using namespace std;
void mqtt()
{
    MQTTScheduler schedule = MQTTScheduler();
    schedule.insert(2, 2, 1);
    schedule.insert(2, 2, 2);
    schedule.insert(4, 1, 3);
    schedule.insert(6, 4, 4);
    cout << "Answer: " << schedule.run() << endl;
    schedule.insert(5, 1, 5);
    cout << "Answer: " << schedule.run() << endl;
    schedule.insert(3, 1, 6);
    cout << "Answer: " << schedule.run() << endl;
    schedule.insert(1, 1, 7);
    cout << "Answer: " << schedule.run() << endl;
    schedule.clear();
    schedule.insert(4, 4, 1);
    schedule.insert(7, 4, 2);
    schedule.insert(6, 2, 3);
    cout << "Answer: " << schedule.run() << endl;
}