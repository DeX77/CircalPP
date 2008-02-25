<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>Alignment.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>Alignment_8cpp</filename>
    <includes id="Alignment_8h" name="Alignment.h" local="yes" imported="no">Alignment.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>Alignment.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>Alignment_8h</filename>
    <namespace>bpp</namespace>
    <namespace>Circal</namespace>
    <class kind="class">Circal::Alignment</class>
  </compound>
  <compound kind="file">
    <name>AlignmentFactory.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>AlignmentFactory_8cpp</filename>
    <includes id="AlignmentFactory_8h" name="AlignmentFactory.h" local="yes" imported="no">AlignmentFactory.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>AlignmentFactory.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>AlignmentFactory_8h</filename>
    <includes id="Alignment_8h" name="Alignment.h" local="yes" imported="no">Alignment.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::AlignmentFactory</class>
    <member kind="typedef">
      <type>std::vector&lt; std::vector&lt; double &gt; &gt;</type>
      <name>ScoreMatrix</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>f94948de4bda3c857d80444f2fbce6ee</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>AlignmentSymbol.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>AlignmentSymbol_8cpp</filename>
    <includes id="AlignmentSymbol_8h" name="AlignmentSymbol.h" local="yes" imported="no">AlignmentSymbol.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>AlignmentSymbol.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>AlignmentSymbol_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::AlignmentSymbol</class>
  </compound>
  <compound kind="file">
    <name>CircalPP.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>CircalPP_8cpp</filename>
    <includes id="MultipleCircularAlignmentFactory_8h" name="MultipleCircularAlignmentFactory.h" local="yes" imported="no">MultipleCircularAlignmentFactory.h</includes>
    <includes id="MultiplePseudoCircularAlignmentFactory_8h" name="MultiplePseudoCircularAlignmentFactory.h" local="yes" imported="no">MultiplePseudoCircularAlignmentFactory.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="CorrectedFasta_8h" name="CorrectedFasta.h" local="yes" imported="no">CorrectedFasta.h</includes>
    <includes id="WhitespaceFasta_8h" name="WhitespaceFasta.h" local="yes" imported="no">WhitespaceFasta.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <includes id="RandomSequence_8h" name="RandomSequence.h" local="yes" imported="no">RandomSequence.h</includes>
    <includes id="VertebrateMitochondrialGenomeAlphabet_8h" name="VertebrateMitochondrialGenomeAlphabet.h" local="yes" imported="no">VertebrateMitochondrialGenomeAlphabet.h</includes>
    <member kind="function">
      <type>std::string</type>
      <name>usage</name>
      <anchorfile>CircalPP_8cpp.html</anchorfile>
      <anchor>10db687c0645469ee3404ed649063237</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>doAllignment</name>
      <anchorfile>CircalPP_8cpp.html</anchorfile>
      <anchor>ea6df4915089872325f1bec671ec55f4</anchor>
      <arglist>(const bpp::Alphabet *alpha, const std::string &amp;seqFilename, const Circal::ScoringModel &amp;scoreM, const std::string &amp;outFilename, const std::string &amp;resultFilename, bool outF, bool resultF, int &amp;delta, bool verbose, bool stepWise)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>CircalPP_8cpp.html</anchorfile>
      <anchor>7be7b3f3b810d259483db57fef9b4c4c</anchor>
      <arglist>(int args, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>CircularAlignmentFactory.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>CircularAlignmentFactory_8cpp</filename>
    <includes id="CircularAlignmentFactory_8h" name="CircularAlignmentFactory.h" local="yes" imported="no">CircularAlignmentFactory.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <includes id="RotatedSequence_8h" name="RotatedSequence.h" local="yes" imported="no">RotatedSequence.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>CircularAlignmentFactory.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>CircularAlignmentFactory_8h</filename>
    <includes id="AlignmentFactory_8h" name="AlignmentFactory.h" local="yes" imported="no">AlignmentFactory.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::CircularAlignmentFactory</class>
  </compound>
  <compound kind="file">
    <name>CorrectedFasta.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>CorrectedFasta_8cpp</filename>
    <includes id="CorrectedFasta_8h" name="CorrectedFasta.h" local="yes" imported="no">CorrectedFasta.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>CorrectedFasta.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>CorrectedFasta_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::CorrectedFasta</class>
  </compound>
  <compound kind="file">
    <name>ErrorClasses.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>ErrorClasses_8cpp</filename>
    <includes id="ErrorClasses_8h" name="ErrorClasses.h" local="yes" imported="no">ErrorClasses.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>ErrorClasses.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>ErrorClasses_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::MatrixOutofBoundsError</class>
  </compound>
  <compound kind="file">
    <name>GenomeAlphabet.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>GenomeAlphabet_8cpp</filename>
    <includes id="GenomeAlphabet_8h" name="GenomeAlphabet.h" local="yes" imported="no">GenomeAlphabet.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>GenomeAlphabet.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>GenomeAlphabet_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::GenomeAlphabet</class>
  </compound>
  <compound kind="file">
    <name>MatrixHelper.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MatrixHelper_8cpp</filename>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="Alignment_8h" name="Alignment.h" local="yes" imported="no">Alignment.h</includes>
    <includes id="PseudoRotatedSequence_8h" name="PseudoRotatedSequence.h" local="yes" imported="no">PseudoRotatedSequence.h</includes>
    <includes id="ErrorClasses_8h" name="ErrorClasses.h" local="yes" imported="no">ErrorClasses.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>MatrixHelper.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MatrixHelper_8h</filename>
    <namespace>bpp</namespace>
    <namespace>Circal</namespace>
    <class kind="class">Circal::MatrixHelper</class>
    <member kind="typedef">
      <type>std::valarray&lt; std::valarray&lt; bool &gt; &gt;</type>
      <name>BoolMatrix</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>b5cbf04ae410080c3875d2691c10d64c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>ScoreMatrix3D</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>8cbdc0f70cb5adc1e0c93eea91ff2883</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>MultipleAlignmentFactory.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MultipleAlignmentFactory_8cpp</filename>
    <includes id="MultipleAlignmentFactory_8h" name="MultipleAlignmentFactory.h" local="yes" imported="no">MultipleAlignmentFactory.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>MultipleAlignmentFactory.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MultipleAlignmentFactory_8h</filename>
    <includes id="AlignmentFactory_8h" name="AlignmentFactory.h" local="yes" imported="no">AlignmentFactory.h</includes>
    <includes id="Alignment_8h" name="Alignment.h" local="yes" imported="no">Alignment.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::MultipleAlignmentFactory</class>
  </compound>
  <compound kind="file">
    <name>MultipleCircularAlignmentFactory.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MultipleCircularAlignmentFactory_8cpp</filename>
    <includes id="MultipleCircularAlignmentFactory_8h" name="MultipleCircularAlignmentFactory.h" local="yes" imported="no">MultipleCircularAlignmentFactory.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>MultipleCircularAlignmentFactory.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MultipleCircularAlignmentFactory_8h</filename>
    <includes id="CircularAlignmentFactory_8h" name="CircularAlignmentFactory.h" local="yes" imported="no">CircularAlignmentFactory.h</includes>
    <includes id="MultipleAlignmentFactory_8h" name="MultipleAlignmentFactory.h" local="yes" imported="no">MultipleAlignmentFactory.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::MultipleCircularAlignmentFactory</class>
  </compound>
  <compound kind="file">
    <name>MultiplePseudoCircularAlignmentFactory.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MultiplePseudoCircularAlignmentFactory_8cpp</filename>
    <includes id="MultiplePseudoCircularAlignmentFactory_8h" name="MultiplePseudoCircularAlignmentFactory.h" local="yes" imported="no">MultiplePseudoCircularAlignmentFactory.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>MultiplePseudoCircularAlignmentFactory.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>MultiplePseudoCircularAlignmentFactory_8h</filename>
    <includes id="PseudoCircularAlignmentFactory_8h" name="PseudoCircularAlignmentFactory.h" local="yes" imported="no">PseudoCircularAlignmentFactory.h</includes>
    <includes id="MultipleAlignmentFactory_8h" name="MultipleAlignmentFactory.h" local="yes" imported="no">MultipleAlignmentFactory.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::MultiplePseudoCircularAlignmentFactory</class>
  </compound>
  <compound kind="file">
    <name>Output.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>Output_8cpp</filename>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <includes id="Alignment_8h" name="Alignment.h" local="yes" imported="no">Alignment.h</includes>
    <includes id="VertebrateMitochondrialGenomeAlphabet_8h" name="VertebrateMitochondrialGenomeAlphabet.h" local="yes" imported="no">VertebrateMitochondrialGenomeAlphabet.h</includes>
    <includes id="PseudoRotatedSequence_8h" name="PseudoRotatedSequence.h" local="yes" imported="no">PseudoRotatedSequence.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>Output.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>Output_8h</filename>
    <namespace>bpp</namespace>
    <namespace>Circal</namespace>
    <class kind="class">Circal::Output</class>
  </compound>
  <compound kind="file">
    <name>PseudoCircularAlignmentFactory.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>PseudoCircularAlignmentFactory_8cpp</filename>
    <includes id="PseudoCircularAlignmentFactory_8h" name="PseudoCircularAlignmentFactory.h" local="yes" imported="no">PseudoCircularAlignmentFactory.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="MatrixHelper_8h" name="MatrixHelper.h" local="yes" imported="no">MatrixHelper.h</includes>
    <includes id="PseudoRotatedSequence_8h" name="PseudoRotatedSequence.h" local="yes" imported="no">PseudoRotatedSequence.h</includes>
    <includes id="Output_8h" name="Output.h" local="yes" imported="no">Output.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>PseudoCircularAlignmentFactory.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>PseudoCircularAlignmentFactory_8h</filename>
    <includes id="AlignmentFactory_8h" name="AlignmentFactory.h" local="yes" imported="no">AlignmentFactory.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::PseudoCircularAlignmentFactory</class>
  </compound>
  <compound kind="file">
    <name>PseudoRotatedSequence.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>PseudoRotatedSequence_8cpp</filename>
    <includes id="PseudoRotatedSequence_8h" name="PseudoRotatedSequence.h" local="yes" imported="no">PseudoRotatedSequence.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>PseudoRotatedSequence.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>PseudoRotatedSequence_8h</filename>
    <includes id="RotatedSequence_8h" name="RotatedSequence.h" local="yes" imported="no">RotatedSequence.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::PseudoRotatedSequence</class>
  </compound>
  <compound kind="file">
    <name>RandomSequence.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>RandomSequence_8cpp</filename>
    <includes id="RandomSequence_8h" name="RandomSequence.h" local="yes" imported="no">RandomSequence.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>RandomSequence.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>RandomSequence_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::RandomSequence</class>
  </compound>
  <compound kind="file">
    <name>RotatedSequence.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>RotatedSequence_8cpp</filename>
    <includes id="RotatedSequence_8h" name="RotatedSequence.h" local="yes" imported="no">RotatedSequence.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>RotatedSequence.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>RotatedSequence_8h</filename>
    <namespace>bpp</namespace>
    <namespace>Circal</namespace>
    <class kind="class">Circal::RotatedSequence</class>
  </compound>
  <compound kind="file">
    <name>ScoringModel.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>ScoringModel_8cpp</filename>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <includes id="AlignmentSymbol_8h" name="AlignmentSymbol.h" local="yes" imported="no">AlignmentSymbol.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>ScoringModel.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>ScoringModel_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::ScoringModel</class>
    <member kind="typedef">
      <type>std::map&lt; const std::string, AlignmentSymbol &gt;</type>
      <name>ModelValues</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>539a9bf369c344e6d500b60b6a109696</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>VertebrateMitochondrialGenomeAlphabet.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>VertebrateMitochondrialGenomeAlphabet_8cpp</filename>
    <includes id="VertebrateMitochondrialGenomeAlphabet_8h" name="VertebrateMitochondrialGenomeAlphabet.h" local="yes" imported="no">VertebrateMitochondrialGenomeAlphabet.h</includes>
    <includes id="AlignmentSymbol_8h" name="AlignmentSymbol.h" local="yes" imported="no">AlignmentSymbol.h</includes>
    <includes id="ScoringModel_8h" name="ScoringModel.h" local="yes" imported="no">ScoringModel.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>VertebrateMitochondrialGenomeAlphabet.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>VertebrateMitochondrialGenomeAlphabet_8h</filename>
    <includes id="GenomeAlphabet_8h" name="GenomeAlphabet.h" local="yes" imported="no">GenomeAlphabet.h</includes>
    <namespace>Circal</namespace>
    <class kind="class">Circal::VertebrateMitochondrialGenomeAlphabet</class>
  </compound>
  <compound kind="file">
    <name>WhitespaceFasta.cpp</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>WhitespaceFasta_8cpp</filename>
    <includes id="WhitespaceFasta_8h" name="WhitespaceFasta.h" local="yes" imported="no">WhitespaceFasta.h</includes>
    <namespace>Circal</namespace>
  </compound>
  <compound kind="file">
    <name>WhitespaceFasta.h</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>WhitespaceFasta_8h</filename>
    <namespace>Circal</namespace>
    <class kind="class">Circal::WhitespaceFasta</class>
  </compound>
  <compound kind="namespace">
    <name>bpp</name>
    <filename>namespacebpp.html</filename>
  </compound>
  <compound kind="namespace">
    <name>Circal</name>
    <filename>namespaceCircal.html</filename>
    <class kind="class">Circal::Alignment</class>
    <class kind="class">Circal::AlignmentFactory</class>
    <class kind="class">Circal::AlignmentSymbol</class>
    <class kind="class">Circal::CircularAlignmentFactory</class>
    <class kind="class">Circal::CorrectedFasta</class>
    <class kind="class">Circal::MatrixOutofBoundsError</class>
    <class kind="class">Circal::GenomeAlphabet</class>
    <class kind="class">Circal::MatrixHelper</class>
    <class kind="class">Circal::MultipleAlignmentFactory</class>
    <class kind="class">Circal::MultipleCircularAlignmentFactory</class>
    <class kind="class">Circal::MultiplePseudoCircularAlignmentFactory</class>
    <class kind="class">Circal::Output</class>
    <class kind="class">Circal::PseudoCircularAlignmentFactory</class>
    <class kind="class">Circal::PseudoRotatedSequence</class>
    <class kind="class">Circal::RandomSequence</class>
    <class kind="class">Circal::RotatedSequence</class>
    <class kind="class">Circal::ScoringModel</class>
    <class kind="class">Circal::VertebrateMitochondrialGenomeAlphabet</class>
    <class kind="class">Circal::WhitespaceFasta</class>
    <member kind="typedef">
      <type>std::vector&lt; std::vector&lt; double &gt; &gt;</type>
      <name>ScoreMatrix</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>f94948de4bda3c857d80444f2fbce6ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::valarray&lt; std::valarray&lt; bool &gt; &gt;</type>
      <name>BoolMatrix</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>b5cbf04ae410080c3875d2691c10d64c</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>ScoreMatrix3D</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>8cbdc0f70cb5adc1e0c93eea91ff2883</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>std::map&lt; const std::string, AlignmentSymbol &gt;</type>
      <name>ModelValues</name>
      <anchorfile>namespaceCircal.html</anchorfile>
      <anchor>539a9bf369c344e6d500b60b6a109696</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::Alignment</name>
    <filename>classCircal_1_1Alignment.html</filename>
    <member kind="function">
      <type></type>
      <name>Alignment</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>6c1d8170276265997eb8729e83a65fa8</anchor>
      <arglist>(const bpp::Alphabet *alpha)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Alignment</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>4a0146b22bceb4dff7a0583378f6145c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>uint</type>
      <name>get_offsetA</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>1ba1a02f0da981faca500e0988bf1ac8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_offsetA</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>222ba2a7ddae54bda89e809b59aa7f12</anchor>
      <arglist>(const uint &amp;orig)</arglist>
    </member>
    <member kind="function">
      <type>uint</type>
      <name>get_offsetB</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>6822c7df4eab6cd256934baa7ae71ee8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_offsetB</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>cdf5f3f56151b4ea4829a826b6b88b63</anchor>
      <arglist>(const uint &amp;orig)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_Score</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>e4dc0a9b6c0b89f3e3c375b963c76ec1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_Score</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>94d66a0aa8ca3a3ecd47801acbbbc16a</anchor>
      <arglist>(const double &amp;s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>addSequence</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>9068711091bc065e281c076da5bd7dcb</anchor>
      <arglist>(const bpp::Sequence *sequence)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>uint</type>
      <name>offsetA</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>8f289b66dbb1373f0d857f5b43f814dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>uint</type>
      <name>offsetB</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>8afb2ec44c7cda3c08759eb544ea822a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>Score</name>
      <anchorfile>classCircal_1_1Alignment.html</anchorfile>
      <anchor>894c310f50b7ec0b44c53d38b03a4c4c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::AlignmentFactory</name>
    <filename>classCircal_1_1AlignmentFactory.html</filename>
    <member kind="function">
      <type></type>
      <name>AlignmentFactory</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>caefc112c9430071c8bcb3761d832d6b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~AlignmentFactory</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>76d219a7899be2c7e3ad545eaec92544</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>NeedlemanWunschAlignment</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>06cf2e0ea959c22e0cd42278b19c94e4</anchor>
      <arglist>(const bpp::Sequence *inA, const bpp::Sequence *inB, const ScoringModel *scoreM, bool verbose=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>GotohAlignment</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>8d6e6e2bb3c061770c10fda190f28877</anchor>
      <arglist>(const bpp::Sequence *inA, const bpp::Sequence *inB, const ScoringModel *scoreM, bool verbose=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>SmithWaterman</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>a16ba19c77b29cd8c2126c8e28e6ee27</anchor>
      <arglist>(const bpp::Sequence *inA, const bpp::Sequence *inB, const ScoringModel *scoreM, bool verbose=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>SmithWatermanAffin</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>fecf05d5484d763107d60794b9510248</anchor>
      <arglist>(const bpp::Sequence *inA, const bpp::Sequence *inB, const ScoringModel *scoreM, bool verbose=false)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>ForwardRecursionNMW</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>23a0c386ae5a5fe2233b21a25bc92bc3</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, ScoreMatrix *D)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingNMW</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>170edf47373ca95d9cd716d4ace79950</anchor>
      <arglist>(const bpp::Sequence *outA, const bpp::Sequence *outB, const ScoringModel *scoreM, const ScoreMatrix *D, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>ForwardRecursionGotoh</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>824141099efd7e6eb36390937a1abec3</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, ScoreMatrix *D, ScoreMatrix *P, ScoreMatrix *Q)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingGotohGlocal</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>a1ac8a1f427263758b520ef65d2d6225</anchor>
      <arglist>(const bpp::Sequence *outA, const bpp::Sequence *outB, const ScoringModel *scoreM, const ScoreMatrix *D, const ScoreMatrix *P, const ScoreMatrix *Q, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingGotohGlobal</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>e9a87f9cf77459708e957d887b2f358d</anchor>
      <arglist>(const bpp::Sequence *outA, const bpp::Sequence *outB, const ScoringModel *scoreM, const ScoreMatrix *D, const ScoreMatrix *P, const ScoreMatrix *Q, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual double</type>
      <name>ForwardRecursionSmithWaterman</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>2554dd57c516fefab4b3c746a10d5a9e</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, ScoreMatrix *D, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingSmithWaterman</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>2c2d383ecd1f50ab572e1813188987e7</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, const ScoreMatrix *D, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual double</type>
      <name>ForwardRecursionSmithWatermanAffin</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>b4cb6cb6ecbffc8d2d6edf46d54b75aa</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, ScoreMatrix *D, ScoreMatrix *P, ScoreMatrix *Q, uint &amp;bi, uint &amp;bj)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingSmithWatermanAffin</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>38e1e5e94033be95ebe5438df2cdae96</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, const ScoreMatrix *D, const ScoreMatrix *P, const ScoreMatrix *Q, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixHelper *</type>
      <name>matrix</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>3686548dbeea38c49ee1cf12c13de560</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Output *</type>
      <name>prettyPrint</name>
      <anchorfile>classCircal_1_1AlignmentFactory.html</anchorfile>
      <anchor>b71d8661bb7da1609fa9cb430204e21c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::AlignmentSymbol</name>
    <filename>classCircal_1_1AlignmentSymbol.html</filename>
    <member kind="function">
      <type></type>
      <name>AlignmentSymbol</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>a5037b9684df598918e945203f517cc1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~AlignmentSymbol</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>4371b3c3d720bc4ed70cdd2b8a564b5a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>symbol</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>c3ca0649bbfcf50c7f50c6ea280a5745</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>type</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>7ed0cee1ba897f9669eb480a453486ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>gapOpen</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>d27fdb207847b400876f66b4b0054211</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>gapExtend</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>7022ea193b1e3f3fba196d9b5026df67</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>match</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>c5e275629a2f83f578f3bab25b1e3682</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>negativeMatch</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>63243049e31d18acafa977a22710e557</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>missmatch</name>
      <anchorfile>classCircal_1_1AlignmentSymbol.html</anchorfile>
      <anchor>b7957127bbfa37693391938563ca601d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::CircularAlignmentFactory</name>
    <filename>classCircal_1_1CircularAlignmentFactory.html</filename>
    <base virtualness="virtual">Circal::AlignmentFactory</base>
    <member kind="function">
      <type></type>
      <name>CircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1CircularAlignmentFactory.html</anchorfile>
      <anchor>ca3d9c621ae46232edae478521d06f6a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1CircularAlignmentFactory.html</anchorfile>
      <anchor>bfcebb9e1c8e590d23943ffed3d9298b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>NeedlemanWunschAlignment</name>
      <anchorfile>classCircal_1_1CircularAlignmentFactory.html</anchorfile>
      <anchor>875288196816426e679eb5bcaf2ec442</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM)</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>GotohAlignment</name>
      <anchorfile>classCircal_1_1CircularAlignmentFactory.html</anchorfile>
      <anchor>2011395ecb3e47fc21987dcc9b9a4812</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::CorrectedFasta</name>
    <filename>classCircal_1_1CorrectedFasta.html</filename>
    <member kind="function">
      <type></type>
      <name>CorrectedFasta</name>
      <anchorfile>classCircal_1_1CorrectedFasta.html</anchorfile>
      <anchor>886f025a1872bf6ed19773f4ecb4406c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CorrectedFasta</name>
      <anchorfile>classCircal_1_1CorrectedFasta.html</anchorfile>
      <anchor>0cf5d4979869682b1217514c7034da39</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>appendFromStream</name>
      <anchorfile>classCircal_1_1CorrectedFasta.html</anchorfile>
      <anchor>00ac2b9dfd12ba34255618acc0193c8c</anchor>
      <arglist>(istream &amp;input, bpp::VectorSequenceContainer &amp;vsc) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::MatrixOutofBoundsError</name>
    <filename>classCircal_1_1MatrixOutofBoundsError.html</filename>
    <member kind="function">
      <type></type>
      <name>MatrixOutofBoundsError</name>
      <anchorfile>classCircal_1_1MatrixOutofBoundsError.html</anchorfile>
      <anchor>11e9285e4f4b82021ccfa5242291a0f3</anchor>
      <arglist>(const uint &amp;start, const uint &amp;end, const uint &amp;size)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~MatrixOutofBoundsError</name>
      <anchorfile>classCircal_1_1MatrixOutofBoundsError.html</anchorfile>
      <anchor>c786252d1a63ba3870ae5c6b36b57346</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::GenomeAlphabet</name>
    <filename>classCircal_1_1GenomeAlphabet.html</filename>
    <member kind="function">
      <type></type>
      <name>GenomeAlphabet</name>
      <anchorfile>classCircal_1_1GenomeAlphabet.html</anchorfile>
      <anchor>1c7ac6e806bb5b5910b34e9af91c9d7a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~GenomeAlphabet</name>
      <anchorfile>classCircal_1_1GenomeAlphabet.html</anchorfile>
      <anchor>60b29009b1ecbbf24af83988d7a28d55</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::MatrixHelper</name>
    <filename>classCircal_1_1MatrixHelper.html</filename>
    <member kind="function">
      <type></type>
      <name>MatrixHelper</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>632fdc53bb0f37a622eba3a4a9cb2632</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~MatrixHelper</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>fb98525415a508eafcf8bc4207b77c13</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ScoreMatrix</type>
      <name>InitializeScoreMatrixDistances</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>69554a977a56437127203bf586e5d594</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMatrix</type>
      <name>InitScoreMatrixWith</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>8c12cd1078ffbb6c4c93ebbb8cb384d5</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const double &amp;init)</arglist>
    </member>
    <member kind="function">
      <type>ScoreMatrix3D</type>
      <name>InitScoreMatrix3DWith</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>d7f73b9f9ada591e7ed652909e240140</anchor>
      <arglist>(const bpp::Sequence *A, const PseudoRotatedSequence *B, const int &amp;delta, const double &amp;init)</arglist>
    </member>
    <member kind="function">
      <type>BoolMatrix</type>
      <name>CreateAdjacenceGraph</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>86ccdcd3fc84854e446918af43d9b9d2</anchor>
      <arglist>(Alignment *pairWiseAlignments, int biggestSequenceSize)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CutRowFromTo</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>baf89eb7d195cd23dbd1b12ddd89ff3d</anchor>
      <arglist>(ScoreMatrix *D, const uint &amp;start, const uint &amp;end)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CutColumnFromTo</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>341c7274b762c4d1c32950816d3c8594</anchor>
      <arglist>(ScoreMatrix *D, const uint &amp;start, const uint &amp;end)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>SearchBestPositionFrom</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>8a261ceaa3de47879c771252a3e38b53</anchor>
      <arglist>(const ScoreMatrix *M, uint &amp;i, uint &amp;j, const ScoringModel *scoreM)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>SearchBestInRow</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>beed35923850dd22a830a22063a77a0a</anchor>
      <arglist>(const ScoreMatrix *M, const ScoringModel *scoreM, const uint &amp;start, const uint &amp;row)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>SearchBestInRow3D</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>0ede39b0287474c056d137f6e812c27f</anchor>
      <arglist>(const ScoreMatrix3D *M, const ScoringModel *scoreM, const uint &amp;start, const uint &amp;row, uint &amp;k)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>SearchBestInColumn</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>e69e461388db01111d662b709d6ab4f2</anchor>
      <arglist>(const ScoreMatrix *M, const ScoringModel *scoreM, const uint &amp;start, const uint &amp;column)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>SearchBestInColumn3D</name>
      <anchorfile>classCircal_1_1MatrixHelper.html</anchorfile>
      <anchor>b3b25325d1f7ee51562ccda3ac015ae6</anchor>
      <arglist>(const ScoreMatrix3D *M, const ScoringModel *scoreM, const uint &amp;start, const uint &amp;column, uint &amp;k)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::MultipleAlignmentFactory</name>
    <filename>classCircal_1_1MultipleAlignmentFactory.html</filename>
    <base virtualness="virtual">Circal::AlignmentFactory</base>
    <member kind="function">
      <type></type>
      <name>MultipleAlignmentFactory</name>
      <anchorfile>classCircal_1_1MultipleAlignmentFactory.html</anchorfile>
      <anchor>377f18b5178ee35a038d734ae9b4f0f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~MultipleAlignmentFactory</name>
      <anchorfile>classCircal_1_1MultipleAlignmentFactory.html</anchorfile>
      <anchor>bfeee62c22ee9a141d09fd45af70ee3c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>GotohalignMultiple</name>
      <anchorfile>classCircal_1_1MultipleAlignmentFactory.html</anchorfile>
      <anchor>7269edee2a627143253ef90711eec760</anchor>
      <arglist>(const bpp::VectorSequenceContainer *input, const ScoringModel *scoreM, bool verbose)</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>NMWalignMultiple</name>
      <anchorfile>classCircal_1_1MultipleAlignmentFactory.html</anchorfile>
      <anchor>12044c37a192e7e9d895d86ec1d426ae</anchor>
      <arglist>(const bpp::VectorSequenceContainer *input, const ScoringModel *scoreM, bool verbose)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::MultipleCircularAlignmentFactory</name>
    <filename>classCircal_1_1MultipleCircularAlignmentFactory.html</filename>
    <base>Circal::CircularAlignmentFactory</base>
    <base>Circal::MultipleAlignmentFactory</base>
    <member kind="function">
      <type></type>
      <name>MultipleCircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1MultipleCircularAlignmentFactory.html</anchorfile>
      <anchor>c2b06487ba86aea2a84ad75d02701590</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~MultipleCircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1MultipleCircularAlignmentFactory.html</anchorfile>
      <anchor>8a67f9c4df17e66d443875693ef320b8</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::MultiplePseudoCircularAlignmentFactory</name>
    <filename>classCircal_1_1MultiplePseudoCircularAlignmentFactory.html</filename>
    <base>Circal::MultipleAlignmentFactory</base>
    <base>Circal::PseudoCircularAlignmentFactory</base>
    <member kind="function">
      <type></type>
      <name>MultiplePseudoCircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1MultiplePseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>35501504926816a4031c14b03b26903b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~MultiplePseudoCircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1MultiplePseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>e4b045519eab95324a91eb9cd2b02d77</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>NMWalignMultiple</name>
      <anchorfile>classCircal_1_1MultiplePseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>47d7db6773b0de0b35387cf0346af2bb</anchor>
      <arglist>(const bpp::VectorSequenceContainer *input, const ScoringModel *scoreM, const int &amp;delta, bool verbose)</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>GotohalignMultiple</name>
      <anchorfile>classCircal_1_1MultiplePseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>e1c5bfaba040571300262dc25c31fb77</anchor>
      <arglist>(const bpp::VectorSequenceContainer *input, const ScoringModel *scoreM, const int &amp;delta, bool verbose)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::Output</name>
    <filename>classCircal_1_1Output.html</filename>
    <member kind="function">
      <type></type>
      <name>Output</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>f7a38f5496d10ae2cfac9646e7c1f6fa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Output</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>e3f937db11ee87443ea1e75049e0706e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>SequencePrettyPrint</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>dea07429ad2ed3d082aa4d9ef85e8114</anchor>
      <arglist>(bpp::Sequence *A)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>AlignmentPrettyPrint</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>ab041c2ce86bb7ab5f682cdab596d06b</anchor>
      <arglist>(Alignment *aln)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>ScoreMatrixPrettyPrint</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>2ed0aa94d088dcffe1c3d1fa633aed0a</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoreMatrix &amp;D)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>AdjacenceMatrixPrettyPrint</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>b92b148d9dcd8b3692140c0d9020ab1b</anchor>
      <arglist>(const BoolMatrix &amp;G)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>TCoffeeLibHeader</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>eeba75375812d32a5ed680a67b406f37</anchor>
      <arglist>(const bpp::VectorSequenceContainer *input)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>TCoffeeAlignFormat</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>8d4d121954ec0c10ad6b3b1268651c08</anchor>
      <arglist>(Alignment *aln, const bpp::VectorSequenceContainer *input)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>TCoffeeLibFooter</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>4f74b8fc9c68d69e1d26ade8c6fbf35d</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>TCoffeeLibFormat</name>
      <anchorfile>classCircal_1_1Output.html</anchorfile>
      <anchor>726f8644afc9acca3f02a3983c7409fe</anchor>
      <arglist>(Alignment *aln, const bpp::VectorSequenceContainer *input)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::PseudoCircularAlignmentFactory</name>
    <filename>classCircal_1_1PseudoCircularAlignmentFactory.html</filename>
    <base virtualness="virtual">Circal::AlignmentFactory</base>
    <member kind="function">
      <type></type>
      <name>PseudoCircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>268a3c8290bfe004aa52a8bb37bcf4d9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~PseudoCircularAlignmentFactory</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>000683a8fc886c342b357bd0cc285854</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>NeedlemanWunschAlignment</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>de2586f4840d130ca24636c2fa3d7c69</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, const int &amp;delta, bool verbose=false)</arglist>
    </member>
    <member kind="function">
      <type>Alignment</type>
      <name>GotohAlignment</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>ebfd497019d9d510398f69824da5c37b</anchor>
      <arglist>(const bpp::Sequence *A, const bpp::Sequence *B, const ScoringModel *scoreM, const int &amp;delta, bool verbose=false)</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual double</type>
      <name>ForwardRecursionSmithWaterman</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>81b0b7b0acdb8ec01af0f9383826a42a</anchor>
      <arglist>(const bpp::Sequence *A, const PseudoRotatedSequence *B, const ScoringModel *scoreM, const int &amp;delta, ScoreMatrix3D *D, uint &amp;bi, uint &amp;bj)</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingSmithWaterman</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>bbebcc77eef75596c7eeb73a24c5d5ed</anchor>
      <arglist>(const bpp::Sequence *A, const PseudoRotatedSequence *B, const ScoringModel *scoreM, const int &amp;delta, const ScoreMatrix3D *D, uint &amp;i, uint &amp;j)</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual double</type>
      <name>ForwardRecursionSmithWatermanAffin</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>74b7860fd17df40cb002981bfa4ceae3</anchor>
      <arglist>(const bpp::Sequence *A, const PseudoRotatedSequence *B, const ScoringModel *scoreM, const int &amp;delta, ScoreMatrix3D *D, ScoreMatrix3D *P, ScoreMatrix3D *Q, uint &amp;bi, uint &amp;bj)</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual Alignment</type>
      <name>BacktrackingSmithWatermanAffin</name>
      <anchorfile>classCircal_1_1PseudoCircularAlignmentFactory.html</anchorfile>
      <anchor>dce787a2350cfd3683dd0682ef064d5a</anchor>
      <arglist>(const bpp::Sequence *A, const PseudoRotatedSequence *B, const ScoringModel *scoreM, const int &amp;delta, const ScoreMatrix3D *D, const ScoreMatrix3D *P, const ScoreMatrix3D *Q, uint &amp;i, uint &amp;j, bool verbose=false)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::PseudoRotatedSequence</name>
    <filename>classCircal_1_1PseudoRotatedSequence.html</filename>
    <base>Circal::RotatedSequence</base>
    <member kind="function">
      <type></type>
      <name>PseudoRotatedSequence</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>bac16ca7b9725950d138748ae91a98b8</anchor>
      <arglist>(const bpp::Sequence *a)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~PseudoRotatedSequence</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>e585ce628e0855bbee4e9ebf5446263f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned int</type>
      <name>size</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>f6a9b351636ffee4468adbf1e58fdb4a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getValue</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>48c964ebdf631807be5e58140e5378c6</anchor>
      <arglist>(unsigned int pos) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual string</type>
      <name>getChar</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>fbfa83464c33febec8196f23cd296a31</anchor>
      <arglist>(unsigned int pos) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const int &amp;</type>
      <name>operator[]</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>56da5d9337b808f82a9e90dcd38f3c31</anchor>
      <arglist>(unsigned int i) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int &amp;</type>
      <name>operator[]</name>
      <anchorfile>classCircal_1_1PseudoRotatedSequence.html</anchorfile>
      <anchor>ae61f22e5caff53b3d0d3b0c2ad2218f</anchor>
      <arglist>(unsigned int i)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::RandomSequence</name>
    <filename>classCircal_1_1RandomSequence.html</filename>
    <member kind="function">
      <type></type>
      <name>RandomSequence</name>
      <anchorfile>classCircal_1_1RandomSequence.html</anchorfile>
      <anchor>41eb65032b007562f49f8cc7acafa8a9</anchor>
      <arglist>(const uint &amp;size, const bpp::Alphabet *alpha)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~RandomSequence</name>
      <anchorfile>classCircal_1_1RandomSequence.html</anchorfile>
      <anchor>cd612e70ea793e398a751d1527244fdb</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::RotatedSequence</name>
    <filename>classCircal_1_1RotatedSequence.html</filename>
    <member kind="function">
      <type></type>
      <name>RotatedSequence</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>40a4c1c343033fd8075f1c77dde8ad68</anchor>
      <arglist>(const bpp::Sequence *a, const uint &amp;i)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~RotatedSequence</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>a7b2c11e48f0c2506e126f8b2a1e3ce9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getValue</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>0b523c27f14abd33c0a5ba92abb9ecfa</anchor>
      <arglist>(unsigned int pos) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual string</type>
      <name>getChar</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>02c66f9b086da552376900b1a2004dcc</anchor>
      <arglist>(unsigned int pos) const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>toString</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>b862eae54b9c8fb912a4ae45a204ab4e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const int &amp;</type>
      <name>operator[]</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>bfa366e9d55c1ffebc5c1dcb2fc11855</anchor>
      <arglist>(unsigned int i) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int &amp;</type>
      <name>operator[]</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>2c55ef48027cac9fa1628d5fca9691f4</anchor>
      <arglist>(unsigned int i)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>uint</type>
      <name>offset</name>
      <anchorfile>classCircal_1_1RotatedSequence.html</anchorfile>
      <anchor>63069c3a23ddd24e52b935ff2fe741d4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::ScoringModel</name>
    <filename>classCircal_1_1ScoringModel.html</filename>
    <member kind="function">
      <type></type>
      <name>ScoringModel</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>6b0d64828fe988c3ac6af739d1382063</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ScoringModel</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>77ac8d2c5ca4db4cfb7beb951524dd3e</anchor>
      <arglist>(const std::string &amp;path)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ScoringModel</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>2dc745eb33b08937bc99c094b4a7230f</anchor>
      <arglist>(std::istream &amp;input)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ScoringModel</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>0f91c4aa719a37596655f30751b485c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>read</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>7d9db7271ea52eeebfb23db5039b0333</anchor>
      <arglist>(const std::string &amp;path)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>read</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>2144c248a25d0e571820eeb4ce727f63</anchor>
      <arglist>(std::istream &amp;input)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>ScoreOf</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>c9093a36aae9269b516c34438d2df279</anchor>
      <arglist>(const std::string &amp;A, const std::string &amp;B) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>ScoreOfGapOpen</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>3f934683024d21df07165eee183ee918</anchor>
      <arglist>(const std::string A) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>ScoreOfGapExtend</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>dc158e70948f9cce9c80ef7df6cdde6f</anchor>
      <arglist>(const std::string A) const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>GetModelSize</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>5a72855024938d033e28072bc5528a1b</anchor>
      <arglist>(void) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>BestOfTwo</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>7a0e45da2b4c5912f08258bd56dd1788</anchor>
      <arglist>(const double A, const double B) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>BestOfThree</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>4a300e1b772adbbda981a5608710c027</anchor>
      <arglist>(const double A, const double B, const double C) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>NormalizeScores</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>9c4029894b6d5d238f1ba870ea59a622</anchor>
      <arglist>(const double &amp;maxValue)</arglist>
    </member>
    <member kind="function">
      <type>ModelValues::const_iterator</type>
      <name>constitStart</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>73e7d5817e567d8731732657be3ae459</anchor>
      <arglist>(void) const </arglist>
    </member>
    <member kind="function">
      <type>ModelValues::const_iterator</type>
      <name>constitEnd</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>f9b1e4968c419cad45035197f800986f</anchor>
      <arglist>(void) const </arglist>
    </member>
    <member kind="function">
      <type>ModelValues::iterator</type>
      <name>itStart</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>5e65c909aadd2a8e559cd36c8169a972</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>ModelValues::iterator</type>
      <name>itEnd</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>544c6c9e67d81f89f554a2db3a8f270f</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ModelValues</type>
      <name>model</name>
      <anchorfile>classCircal_1_1ScoringModel.html</anchorfile>
      <anchor>ede10aaf341f9809bb34108f08977ce0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::VertebrateMitochondrialGenomeAlphabet</name>
    <filename>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</filename>
    <base>Circal::GenomeAlphabet</base>
    <member kind="function">
      <type></type>
      <name>VertebrateMitochondrialGenomeAlphabet</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>151482f92385b56201bb05ea9b583a28</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>VertebrateMitochondrialGenomeAlphabet</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>4df3659abd819e42f5b9772564b1a0bb</anchor>
      <arglist>(const ScoringModel *scm)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~VertebrateMitochondrialGenomeAlphabet</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>8a3ec2256f7e2245aa296a185b801015</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getSize</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>8fbac4c626c774463bca276315584849</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getNumberOfTypes</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>6183076c0da262259573e3f2c200851f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getUnknownCharacterCode</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>3301033c95c2504db33946db676d93f8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>getAlphabetType</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>475afe4d4419ae921934aa8b70ff8b59</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isUnresolved</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>790b4cbe060af422360169015f33c91f</anchor>
      <arglist>(int state) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isUnresolved</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>8633d1fb81913be79c8fdb439be574d3</anchor>
      <arglist>(const string &amp;state) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convertModeltoAlphabet</name>
      <anchorfile>classCircal_1_1VertebrateMitochondrialGenomeAlphabet.html</anchorfile>
      <anchor>6b3fa1f62294d8c683be6e705f963bce</anchor>
      <arglist>(const ScoringModel *scm)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Circal::WhitespaceFasta</name>
    <filename>classCircal_1_1WhitespaceFasta.html</filename>
    <member kind="function">
      <type></type>
      <name>WhitespaceFasta</name>
      <anchorfile>classCircal_1_1WhitespaceFasta.html</anchorfile>
      <anchor>1f47805c3e0292dc6b7be224a3ec06bf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~WhitespaceFasta</name>
      <anchorfile>classCircal_1_1WhitespaceFasta.html</anchorfile>
      <anchor>9cde8d78c3f94f7e599fed3dd4aab2d1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>appendFromStream</name>
      <anchorfile>classCircal_1_1WhitespaceFasta.html</anchorfile>
      <anchor>8e35d785e2d97a71b081bb7a02c12309</anchor>
      <arglist>(istream &amp;input, bpp::VectorSequenceContainer &amp;vsc) const </arglist>
    </member>
  </compound>
  <compound kind="dir">
    <name>Uni/workspace/CircalPP/</name>
    <path>/home/dex/Uni/workspace/CircalPP/</path>
    <filename>dir_a1efb2ef5cbe68894b97b30da26ac79f.html</filename>
    <dir>Uni/workspace/CircalPP/src/</dir>
  </compound>
  <compound kind="dir">
    <name>Uni/workspace/CircalPP/src/</name>
    <path>/home/dex/Uni/workspace/CircalPP/src/</path>
    <filename>dir_bce082e1b219ca11578ca2e11a91790f.html</filename>
    <file>Alignment.cpp</file>
    <file>Alignment.h</file>
    <file>AlignmentFactory.cpp</file>
    <file>AlignmentFactory.h</file>
    <file>AlignmentSymbol.cpp</file>
    <file>AlignmentSymbol.h</file>
    <file>CircalPP.cpp</file>
    <file>CircularAlignmentFactory.cpp</file>
    <file>CircularAlignmentFactory.h</file>
    <file>CorrectedFasta.cpp</file>
    <file>CorrectedFasta.h</file>
    <file>ErrorClasses.cpp</file>
    <file>ErrorClasses.h</file>
    <file>GenomeAlphabet.cpp</file>
    <file>GenomeAlphabet.h</file>
    <file>MatrixHelper.cpp</file>
    <file>MatrixHelper.h</file>
    <file>MultipleAlignmentFactory.cpp</file>
    <file>MultipleAlignmentFactory.h</file>
    <file>MultipleCircularAlignmentFactory.cpp</file>
    <file>MultipleCircularAlignmentFactory.h</file>
    <file>MultiplePseudoCircularAlignmentFactory.cpp</file>
    <file>MultiplePseudoCircularAlignmentFactory.h</file>
    <file>Output.cpp</file>
    <file>Output.h</file>
    <file>PseudoCircularAlignmentFactory.cpp</file>
    <file>PseudoCircularAlignmentFactory.h</file>
    <file>PseudoRotatedSequence.cpp</file>
    <file>PseudoRotatedSequence.h</file>
    <file>RandomSequence.cpp</file>
    <file>RandomSequence.h</file>
    <file>RotatedSequence.cpp</file>
    <file>RotatedSequence.h</file>
    <file>ScoringModel.cpp</file>
    <file>ScoringModel.h</file>
    <file>VertebrateMitochondrialGenomeAlphabet.cpp</file>
    <file>VertebrateMitochondrialGenomeAlphabet.h</file>
    <file>WhitespaceFasta.cpp</file>
    <file>WhitespaceFasta.h</file>
  </compound>
  <compound kind="dir">
    <name>Uni/</name>
    <path>/home/dex/Uni/</path>
    <filename>dir_116dda80bdb57e721ab4ea70d7e41a47.html</filename>
    <dir>Uni/workspace/</dir>
  </compound>
  <compound kind="dir">
    <name>Uni/workspace/</name>
    <path>/home/dex/Uni/workspace/</path>
    <filename>dir_b6ef3898a1f2d8a7de5f0f56f9b40465.html</filename>
    <dir>Uni/workspace/CircalPP/</dir>
  </compound>
</tagfile>
