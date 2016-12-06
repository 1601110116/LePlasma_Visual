/*
 * ExternalField.h
 *
 *  Created on: Mar 13, 2016
 *      Author: zlstudio
 */

#ifndef EXTERNALFIELD_H_
#define EXTERNALFIELD_H_

#include "Grid.h"

class ExternalField {
	public:
		ExternalField();
		virtual ~ExternalField();

		static void setExternalField(Grid* grid);
};

#endif /* EXTERNALFIELD_H_ */
