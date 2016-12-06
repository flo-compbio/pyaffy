..
    Copyright (c) 2016 Florian Wagner
    
    This file is part of pyAffy.
    
    pyAffy is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License, Version 3,
    as published by the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

Changelog
=========

Version 0.3.0 (2016-03-05)
--------------------------

- Improved estimation of RMA background parameter mu (mean background value).
  The previous procedure was too inflexible to handle arbitrary intensity
  ranges.

- Added a parser for the "Command Console" CEL file format.

- CEL file parsers now use ISO-8859-1 encoding for text values.


Version 0.3.2 (2016-12-06)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- Added Python 3 compatibility