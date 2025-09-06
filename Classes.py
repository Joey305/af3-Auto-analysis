class Atom:

	def __init__(self,pdbfileline='', twoCharacterChain=0, legalize=0):

		self.pdbfileline = pdbfileline
		self.record_name = ""
		self.atom_serial, self.hex_atom_serial = 0, 0
		self.atom_name = ""
		self.alternate_location = ""
		self.residue_name = ""
		self.chain_identifier = ""
		self.residue_sequence_number, self.hex_residue_sequence_number = 0, 0
		self.code_for_insertion_residues = ""
		self.x = 0.
		self.y = 0.
		self.z = 0.
		self.occupancy = 0.
		self.temperature_factor = 0.
		self.segment_identifier = ""
		self.element_symbol = ""
		self.charge = ""
		self.atom_key = ""
		self.residue_key = ""
		self.pdbcode = ""

		if self.pdbfileline:
			
			line = self.pdbfileline.strip()
			self.record_name = line[0:6]  
			try:
				# < 99,999 atoms in the PDB file.
				#################################
				self.atom_serial = int(line[6:11]) 
			except:
				# > 99,999 atoms in the PDB file.
				# After 99,999 atoms numbering is hexidecimal.
				##############################################
				self.atom_serial = int(line[6:11], 16)
				self.hex_atom_serial = 1

			self.atom_name = line[12:16].split()[0]
			
			# Hydrogens need four character spaces.
			# All other atoms only need at most three character spaces.
			# In the case of non-hydrogen atoms, add a space preceding the atom name.
			# This preserves the proper PDB format when and if the Atom() instance is written to file.
			##########################################################################################
			if len(self.atom_name) != 4:
				self.atom_name = " " + self.atom_name
					
			self.alternate_location = line[16:17]
			if self.alternate_location == " ":
				self.alternate_location = ""
			self.residue_name = line[17:20] 
			#Remove any white-space in the name of the residue.
			temp_residue_name = ""
			for char in self.residue_name:
				if char != " ":
					temp_residue_name += char
			self.residue_name = temp_residue_name 

			# Accept 1 and 2-character chains by argument.
			##############################################
			if twoCharacterChain:
				# 2-character chain.
				#############################
				self.chain_identifier = line[21:23]
			else:
				# 1-character chain.
				##########################
				self.chain_identifier = line[21:22]

			if self.chain_identifier == " ":
				self.chain_identifier = "NULL"
			else:
				# Remove whitespace if necessary.
				#################################
				self.chain_identifier = self.chain_identifier.split()[0]

			# Accepts 3 and 4-character residue numbers by argument.
			########################################################
			if twoCharacterChain:
				try:
					# int() with base 10.
					#####################
					self.residue_sequence_number = int(line[23:26])
				except:
					# residue is numbered using a hexidecimal system.
					#################################################
					self.residue_sequence_number = int(line[23:26], 16)
					self.hex_residue_sequence_number = 1
			else:
				try:
					self.residue_sequence_number = int(line[22:26])
				except:
					# > 999 residues in the PDB file.
					# After 999 residues numbering is hexidecimal.
					##############################################
					self.residue_sequence_number = int(line[22:26], 16)
					self.hex_residue_sequence_number = 1

			self.code_for_insertion_residues = line[26:27]
			self.x = float(line[30:38]) 
			self.y = float(line[38:46])
			self.z = float(line[46:54])
			if line[54:60]:
				self.occupancy = float(line[54:60])
			if line[60:66]:
				self.temperature_factor = float(line[60:66])
			# XXXX 2019.04.10 shifted the next three lines left...test
			self.segment_identifier = line[72:76] 
			self.element_symbol = line[76:78]
			self.charge = line[78:79]

			self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
			self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)

	def __repr__(self, forceHex=0):
				
		###############################################
		# If the conditional argument forceHex is true,
		# the string respresentation of the atom serial 
		# number and residue number will be written in 
		# hexidecimal format.
		###############################################

		if self.chain_identifier == "NULL":
			chain_identifier = " "
		else:
			chain_identifier = self.chain_identifier
				
		# Robust representation of atom_serial.
		#######################################
		stringAtomSerial = ""
		if self.hex_atom_serial and forceHex:
			stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
		else:
			stringAtomSerial = "%5i" % self.atom_serial

		# Robust representation of residue_sequence_number.
		####################################################   
		stringResidueNumber = ""
		if self.hex_residue_sequence_number and forceHex:
			stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
		else:
			stringResidueNumber = "%3i" % self.residue_sequence_number

		# Format the atom data in standard PDB file format.
		###################################################
		pdbformat  = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
		# Check for 2-character chain identifier.
		# If so, adjust chain and residue number fields accordingly.
		############################################################
		if len(chain_identifier) > 1:
			pdbformat  = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
		pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n"

		return pdbformat % (self.record_name, 
							stringAtomSerial, 
							self.atom_name,
							self.alternate_location,
							self.residue_name,
							chain_identifier,
							stringResidueNumber,
							self.code_for_insertion_residues,
							self.x,
							self.y,
							self.z,
							self.occupancy,
							self.temperature_factor,
							self.segment_identifier,
							self.element_symbol,
							self.charge)

class PseudoAtom:

	def __init__(self):

		self.record_name = "ATOM  "
		self.atom_serial, self.hex_atom_serial = 1, 0
		self.atom_name = " O  " #one space to the left, two to the right
		self.alternate_location = ""
		self.residue_name = "PSA"
		self.chain_identifier = "A"
		self.residue_sequence_number, self.hex_residue_sequence_number = 1, 0
		self.code_for_insertion_residues = ""
		self.x = 0.
		self.y = 0.
		self.z = 0.
		self.occupancy = 0.
		self.temperature_factor = 0.
		self.segment_identifier = ""
		self.element_symbol = ""
		self.charge = ""
		self.atom_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier, self.atom_name)
		self.residue_key = (self.residue_sequence_number, self.residue_name, self.chain_identifier)
		self.pdbcode = ""

	def __repr__(self, forceHex=0):
				
		###############################################
		# If the conditional argument forceHex is true,
		# the string respresentation of the atom serial 
		# number and residue number will be written in 
		# hexidecimal format.
		###############################################

		if self.chain_identifier == "NULL":
			chain_identifier = " "
		else:
			chain_identifier = self.chain_identifier
				
		# Robust representation of atom_serial.
		#######################################
		stringAtomSerial = ""
		if self.hex_atom_serial and forceHex:
			stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
		else:
			stringAtomSerial = "%5i" % self.atom_serial

		# Robust representation of residue_sequence_number.
		####################################################   
		stringResidueNumber = ""
		if self.hex_residue_sequence_number and forceHex:
			stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
		else:
			stringResidueNumber = "%3i" % self.residue_sequence_number

		# Format the atom data in standard PDB file format.
		###################################################
		pdbformat  = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
		# Check for 2-character chain identifier.
		# If so, adjust chain and residue number fields accordingly.
		############################################################
		if len(chain_identifier) > 1:
			pdbformat  = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
		pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n"

		return pdbformat % (self.record_name, 
							stringAtomSerial, 
							self.atom_name,
							self.alternate_location,
							self.residue_name,
							chain_identifier,
							stringResidueNumber,
							self.code_for_insertion_residues,
							self.x,
							self.y,
							self.z,
							self.occupancy,
							self.temperature_factor,
							self.segment_identifier,
							self.element_symbol,
							self.charge)