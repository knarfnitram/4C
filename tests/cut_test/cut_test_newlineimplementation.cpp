// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_combintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

// Tests the new line implementation
//(to sides with 3 common cut points - because of tolerances!)
void test_christoph_1()
{
  Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.12128896315179579779020002661127;
    tri3_xyze(1, 0) = 1.2847875528385353201941032792732e-12;
    tri3_xyze(2, 0) = 0.84523634915106971021714343805797;
    nids.push_back(18219);
    tri3_xyze(0, 1) = -0.1267913382068394101409580798645;
    tri3_xyze(1, 1) = 1.0054029258101650565063525655241e-12;
    tri3_xyze(2, 1) = 0.84221139326627747490050523992977;
    nids.push_back(18215);
    tri3_xyze(0, 2) = -0.12192686786859271985683506045461;
    tri3_xyze(1, 2) = -0.016051976670389100920743885581032;
    tri3_xyze(2, 2) = 0.84372387122994985109158960767672;
    nids.push_back(-103);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.1267913382068394101409580798645;
    tri3_xyze(1, 0) = 1.0054029258101650565063525655241e-12;
    tri3_xyze(2, 0) = 0.84221139326627747490050523992977;
    nids.push_back(18215);
    tri3_xyze(0, 1) = -0.12247102814923389868528857959973;
    tri3_xyze(1, 1) = -0.032816013071493593811212008404254;
    tri3_xyze(2, 1) = 0.84221139330966088287766524445033;
    nids.push_back(18213);
    tri3_xyze(0, 2) = -0.12192686786859271985683506045461;
    tri3_xyze(1, 2) = -0.016051976670389100920743885581032;
    tri3_xyze(2, 2) = 0.84372387122994985109158960767672;
    nids.push_back(-103);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.12247102814923389868528857959973;
    tri3_xyze(1, 0) = -0.032816013071493593811212008404254;
    tri3_xyze(2, 0) = 0.84221139330966088287766524445033;
    nids.push_back(18213);
    tri3_xyze(0, 1) = -0.1267913382068394101409580798645;
    tri3_xyze(1, 1) = 1.0054029258101650565063525655241e-12;
    tri3_xyze(2, 1) = 0.84221139326627747490050523992977;
    nids.push_back(18215);
    tri3_xyze(0, 2) = -0.12712784003924842979316167657089;
    tri3_xyze(1, 2) = -0.016736697646654755816664206236055;
    tri3_xyze(2, 2) = 0.8403660264535651736039767456532;
    nids.push_back(-100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13187119793824150781880177873973;
    tri3_xyze(1, 0) = 6.8381858837254690594033811057439e-13;
    tri3_xyze(2, 0) = 0.83852065959948562934300753113348;
    nids.push_back(18211);
    tri3_xyze(0, 1) = -0.12737779586267888864981046026514;
    tri3_xyze(1, 1) = -0.034130777516814647665199800030678;
    tri3_xyze(2, 1) = 0.83852065963883659627242650458356;
    nids.push_back(18209);
    tri3_xyze(0, 2) = -0.12712784003924842979316167657089;
    tri3_xyze(1, 2) = -0.016736697646654755816664206236055;
    tri3_xyze(2, 2) = 0.8403660264535651736039767456532;
    nids.push_back(-100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13566510187142094556733695753792;
    tri3_xyze(1, 0) = 0.0363513544683356296105536387131;
    tri3_xyze(2, 0) = 0.82938426095883910349471079825889;
    nids.push_back(21391);
    tri3_xyze(0, 1) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 1) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 1) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 2) = -0.13609086066689141114594008286076;
    tri3_xyze(1, 2) = 0.017916701699260025038018540044504;
    tri3_xyze(2, 2) = 0.83180330707877470075573000940494;
    nids.push_back(-96);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13644842991250552555193564785441;
    tri3_xyze(1, 0) = 3.0805900593529961445978665457446e-13;
    tri3_xyze(2, 0) = 0.83422235319976034695343969360692;
    nids.push_back(18207);
    tri3_xyze(0, 1) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 1) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 1) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 2) = -0.13609086066118811220171380682586;
    tri3_xyze(1, 2) = -0.01791670169757834940593710371104;
    tri3_xyze(2, 2) = 0.83180330707815763879864334739977;
    nids.push_back(-94);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 0) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 0) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 1) = -0.13644842991250552555193564785441;
    tri3_xyze(1, 1) = 3.0805900593529961445978665457446e-13;
    tri3_xyze(2, 1) = 0.83422235319976034695343969360692;
    nids.push_back(18207);
    tri3_xyze(0, 2) = -0.13609086066689141114594008286076;
    tri3_xyze(1, 2) = 0.017916701699260025038018540044504;
    tri3_xyze(2, 2) = 0.83180330707877470075573000940494;
    nids.push_back(-96);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 0) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 0) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 1) = -0.13566510187142094556733695753792;
    tri3_xyze(1, 1) = 0.0363513544683356296105536387131;
    tri3_xyze(2, 1) = 0.82938426095883910349471079825889;
    nids.push_back(21391);
    tri3_xyze(0, 2) = -0.13971155692836079165175533489673;
    tri3_xyze(1, 2) = 0.018393375405183418069832157470955;
    tri3_xyze(2, 2) = 0.82673347166320310108744706667494;
    nids.push_back(-93);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13187119793824150781880177873973;
    tri3_xyze(1, 0) = 6.8381858837254690594033811057439e-13;
    tri3_xyze(2, 0) = 0.83852065959948562934300753113348;
    nids.push_back(18211);
    tri3_xyze(0, 1) = -0.13644842991250552555193564785441;
    tri3_xyze(1, 1) = 3.0805900593529961445978665457446e-13;
    tri3_xyze(2, 1) = 0.83422235319976034695343969360692;
    nids.push_back(18207);
    tri3_xyze(0, 2) = -0.13187412153870492481644305371447;
    tri3_xyze(1, 2) = -0.017361557460313127576601033297266;
    tri3_xyze(2, 2) = 0.83637150641930801953094487544149;
    nids.push_back(-97);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 0) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 0) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 1) = -0.13566510185868788096463788406254;
    tri3_xyze(1, 1) = -0.036351354465087172551651661933647;
    tri3_xyze(2, 1) = 0.82938426095735451326618203893304;
    nids.push_back(18201);
    tri3_xyze(0, 2) = -0.13609086066118811220171380682586;
    tri3_xyze(1, 2) = -0.01791670169757834940593710371104;
    tri3_xyze(2, 2) = 0.83180330707815763879864334739977;
    nids.push_back(-94);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13891494439376805836161565821385;
    tri3_xyze(1, 0) = -0.037222147149526534082308870665656;
    tri3_xyze(2, 0) = 0.82408268240795645720453421745333;
    nids.push_back(18197);
    tri3_xyze(0, 1) = -0.13566510185868788096463788406254;
    tri3_xyze(1, 1) = -0.036351354465087172551651661933647;
    tri3_xyze(2, 1) = 0.82938426095735451326618203893304;
    nids.push_back(18201);
    tri3_xyze(0, 2) = -0.13971155692116099533706119473209;
    tri3_xyze(1, 2) = -0.018393375403819994806653426167031;
    tri3_xyze(2, 2) = 0.82673347166241861749824693106348;
    nids.push_back(-91);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13566510185868788096463788406254;
    tri3_xyze(1, 0) = -0.036351354465087172551651661933647;
    tri3_xyze(2, 0) = 0.82938426095735451326618203893304;
    nids.push_back(18201);
    tri3_xyze(0, 1) = -0.13179906244139383275637555925641;
    tri3_xyze(1, 1) = -0.035315452325429728952510544104371;
    tri3_xyze(2, 1) = 0.83422235323914961657720823495765;
    nids.push_back(18205);
    tri3_xyze(0, 2) = -0.13609086066118811220171380682586;
    tri3_xyze(1, 2) = -0.01791670169757834940593710371104;
    tri3_xyze(2, 2) = 0.83180330707815763879864334739977;
    nids.push_back(-94);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 0) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 0) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 1) = -0.14381533300002277697693386926403;
    tri3_xyze(1, 1) = -5.6172194965979902003090033085764e-13;
    tri3_xyze(2, 1) = 0.82408268236799719907992312073475;
    nids.push_back(18199);
    tri3_xyze(0, 2) = -0.13971155692116099533706119473209;
    tri3_xyze(1, 2) = -0.018393375403819994806653426167031;
    tri3_xyze(2, 2) = 0.82673347166241861749824693106348;
    nids.push_back(-91);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.12737779586267888864981046026514;
    tri3_xyze(1, 0) = -0.034130777516814647665199800030678;
    tri3_xyze(2, 0) = 0.83852065963883659627242650458356;
    nids.push_back(18209);
    tri3_xyze(0, 1) = -0.12247102814923389868528857959973;
    tri3_xyze(1, 1) = -0.032816013071493593811212008404254;
    tri3_xyze(2, 1) = 0.84221139330966088287766524445033;
    nids.push_back(18213);
    tri3_xyze(0, 2) = -0.12712784003924842979316167657089;
    tri3_xyze(1, 2) = -0.016736697646654755816664206236055;
    tri3_xyze(2, 2) = 0.8403660264535651736039767456532;
    nids.push_back(-100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.14381533300002277697693386926403;
    tri3_xyze(1, 0) = -5.6172194965979902003090033085764e-13;
    tri3_xyze(2, 0) = 0.82408268236799719907992312073475;
    nids.push_back(18199);
    tri3_xyze(0, 1) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 1) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 1) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 2) = -0.13971155692836079165175533489673;
    tri3_xyze(1, 2) = 0.018393375405183418069832157470955;
    tri3_xyze(2, 2) = 0.82673347166320310108744706667494;
    nids.push_back(-93);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.12737779586267888864981046026514;
    tri3_xyze(1, 0) = -0.034130777516814647665199800030678;
    tri3_xyze(2, 0) = 0.83852065963883659627242650458356;
    nids.push_back(18209);
    tri3_xyze(0, 1) = -0.13187119793824150781880177873973;
    tri3_xyze(1, 1) = 6.8381858837254690594033811057439e-13;
    tri3_xyze(2, 1) = 0.83852065959948562934300753113348;
    nids.push_back(18211);
    tri3_xyze(0, 2) = -0.13187412153870492481644305371447;
    tri3_xyze(1, 2) = -0.017361557460313127576601033297266;
    tri3_xyze(2, 2) = 0.83637150641930801953094487544149;
    nids.push_back(-97);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.14381533300002277697693386926403;
    tri3_xyze(1, 0) = -5.6172194965979902003090033085764e-13;
    tri3_xyze(2, 0) = 0.82408268236799719907992312073475;
    nids.push_back(18199);
    tri3_xyze(0, 1) = -0.13891494439376805836161565821385;
    tri3_xyze(1, 1) = -0.037222147149526534082308870665656;
    tri3_xyze(2, 1) = 0.82408268240795645720453421745333;
    nids.push_back(18197);
    tri3_xyze(0, 2) = -0.13971155692116099533706119473209;
    tri3_xyze(1, 2) = -0.018393375403819994806653426167031;
    tri3_xyze(2, 2) = 0.82673347166241861749824693106348;
    nids.push_back(-91);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13179906244139383275637555925641;
    tri3_xyze(1, 0) = -0.035315452325429728952510544104371;
    tri3_xyze(2, 0) = 0.83422235323914961657720823495765;
    nids.push_back(18205);
    tri3_xyze(0, 1) = -0.13644842991250552555193564785441;
    tri3_xyze(1, 1) = 3.0805900593529961445978665457446e-13;
    tri3_xyze(2, 1) = 0.83422235319976034695343969360692;
    nids.push_back(18207);
    tri3_xyze(0, 2) = -0.13609086066118811220171380682586;
    tri3_xyze(1, 2) = -0.01791670169757834940593710371104;
    tri3_xyze(2, 2) = 0.83180330707815763879864334739977;
    nids.push_back(-94);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13179906244139383275637555925641;
    tri3_xyze(1, 0) = -0.035315452325429728952510544104371;
    tri3_xyze(2, 0) = 0.83422235323914961657720823495765;
    nids.push_back(18205);
    tri3_xyze(0, 1) = -0.12737779586267888864981046026514;
    tri3_xyze(1, 1) = -0.034130777516814647665199800030678;
    tri3_xyze(2, 1) = 0.83852065963883659627242650458356;
    nids.push_back(18209);
    tri3_xyze(0, 2) = -0.13187412153870492481644305371447;
    tri3_xyze(1, 2) = -0.017361557460313127576601033297266;
    tri3_xyze(2, 2) = 0.83637150641930801953094487544149;
    nids.push_back(-97);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13566510185868788096463788406254;
    tri3_xyze(1, 0) = -0.036351354465087172551651661933647;
    tri3_xyze(2, 0) = 0.82938426095735451326618203893304;
    nids.push_back(18201);
    tri3_xyze(0, 1) = -0.14045084843216529280063298301684;
    tri3_xyze(1, 1) = -1.045547837172444238882512701494e-13;
    tri3_xyze(2, 1) = 0.82938426091636585635313849707018;
    nids.push_back(18203);
    tri3_xyze(0, 2) = -0.13971155692116099533706119473209;
    tri3_xyze(1, 2) = -0.018393375403819994806653426167031;
    tri3_xyze(2, 2) = 0.82673347166241861749824693106348;
    nids.push_back(-91);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.13644842991250552555193564785441;
    tri3_xyze(1, 0) = 3.0805900593529961445978665457446e-13;
    tri3_xyze(2, 0) = 0.83422235319976034695343969360692;
    nids.push_back(18207);
    tri3_xyze(0, 1) = -0.13179906244139383275637555925641;
    tri3_xyze(1, 1) = -0.035315452325429728952510544104371;
    tri3_xyze(2, 1) = 0.83422235323914961657720823495765;
    nids.push_back(18205);
    tri3_xyze(0, 2) = -0.13187412153870492481644305371447;
    tri3_xyze(1, 2) = -0.017361557460313127576601033297266;
    tri3_xyze(2, 2) = 0.83637150641930801953094487544149;
    nids.push_back(-97);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.1267913382068394101409580798645;
    tri3_xyze(1, 0) = 1.0054029258101650565063525655241e-12;
    tri3_xyze(2, 0) = 0.84221139326627747490050523992977;
    nids.push_back(18215);
    tri3_xyze(0, 1) = -0.13187119793824150781880177873973;
    tri3_xyze(1, 1) = 6.8381858837254690594033811057439e-13;
    tri3_xyze(2, 1) = 0.83852065959948562934300753113348;
    nids.push_back(18211);
    tri3_xyze(0, 2) = -0.12712784003924842979316167657089;
    tri3_xyze(1, 2) = -0.016736697646654755816664206236055;
    tri3_xyze(2, 2) = 0.8403660264535651736039767456532;
    nids.push_back(-100);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.125;
    hex8_xyze(1, 0) = -0.025000000000000022204460492503131;
    hex8_xyze(2, 0) = 0.82499999999999984456877655247808;
    nids.push_back(14877);
    hex8_xyze(0, 1) = -0.12500000000000002775557561562891;
    hex8_xyze(1, 1) = 0;
    hex8_xyze(2, 1) = 0.82499999999999984456877655247808;
    nids.push_back(14878);
    hex8_xyze(0, 2) = -0.1499999999999999666933092612453;
    hex8_xyze(1, 2) = 0;
    hex8_xyze(2, 2) = 0.82499999999999984456877655247808;
    nids.push_back(14899);
    hex8_xyze(0, 3) = -0.15000000000000002220446049250313;
    hex8_xyze(1, 3) = -0.025000000000000022204460492503131;
    hex8_xyze(2, 3) = 0.82499999999999984456877655247808;
    nids.push_back(14898);
    hex8_xyze(0, 4) = -0.125;
    hex8_xyze(1, 4) = -0.025000000000000008326672684688674;
    hex8_xyze(2, 4) = 0.85000000000000008881784197001252;
    nids.push_back(15318);
    hex8_xyze(0, 5) = -0.12500000000000002775557561562891;
    hex8_xyze(1, 5) = 6.9388939039072283776476979255676e-18;
    hex8_xyze(2, 5) = 0.84999999999999997779553950749687;
    nids.push_back(15319);
    hex8_xyze(0, 6) = -0.1499999999999999666933092612453;
    hex8_xyze(1, 6) = 6.9388939039072283776476979255676e-18;
    hex8_xyze(2, 6) = 0.84999999999999986677323704498122;
    nids.push_back(15340);
    hex8_xyze(0, 7) = -0.15000000000000002220446049250313;
    hex8_xyze(1, 7) = -0.025000000000000008326672684688674;
    hex8_xyze(2, 7) = 0.84999999999999986677323704498122;
    nids.push_back(15339);

    intersection.add_element(1, nids, hex8_xyze, Core::FE::CellType::hex8);

    intersection.cut_test_cut(true, Cut::VCellGaussPts_DirectDivergence);
  }

  std::vector<double> dirDivVol;

  Cut::Mesh mesh = intersection.normal_mesh();
  const std::list<std::shared_ptr<Cut::VolumeCell>>& other_cells = mesh.volume_cells();
  for (std::list<std::shared_ptr<Cut::VolumeCell>>::const_iterator i = other_cells.begin();
      i != other_cells.end(); ++i)
  {
    Cut::VolumeCell* vc = &**i;
    std::cout << "Volume of Volumecell: " << vc->volume() << std::endl;
  }
}
