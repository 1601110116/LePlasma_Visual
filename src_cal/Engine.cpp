/*
 * Engine.cpp
 *
 *  Created on: 2016-5-19
 *      Author: zlstudio
 */

#include <Engine.h>

Engine::Engine(Grid* _grid,double dt):grid(_grid),deltaT(dt) {
}

Engine::~Engine() {
}

void Engine::update(const Range&){}
