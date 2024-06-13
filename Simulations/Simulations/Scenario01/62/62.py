#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.90317181188022]),
array([29.7452689244122, 25.009508901748152, 24.30639715010624, 23.900580810013714, 27.740168568136568, 21.240386023078667, 22.240994166180215, 20.860424073022898, 22.16015051115806, 21.827300514785108, 19.76562318106745, 19.99672765384423]),
array([31.089469996370187, 23.94012896092908, 27.219396699900972, 23.8378503995662, 24.705429714316292, 24.236294954663624, 27.12721359296346, 24.091594226895054, 26.373474547738258, 25.457584902408193, 20.53060303860059, 20.864537720123234, 22.07417181323089, 21.611683864518398, 21.870497531861734, 21.300484362963644, 21.57289376314075, 17.51227913979323, 18.385007170719415, 20.402533392823077, 20.41871867577053, 15.309134748616408, 13.936367864406817, 15.165962368386767, 11.872654862091917, 9.909995322312152, 11.191918715080984, 9.651486395716873, 11.208353463035731, 9.812889792042132, 9.584597037708606, 9.030005544743025, 11.226409989616828]),
array([29.612787972483485, 28.44697450636056, 28.919561572398845, 24.14551699955641, 27.237577335134063, 27.23844484333002, 24.683259015460834, 24.821385844328564, 27.936417714442914, 25.021486711509255, 26.59984474934182]),
array([28.12918709964024, 28.72492150261483, 28.808713570930433, 28.625606032135014, 29.032194110611336, 23.03974598044694, 23.87651927926251, 27.285937361247175, 23.802901730959505, 23.28652343284302, 24.528873147122418, 26.412428155416, 27.66812047955604, 25.871828038482228, 26.065122107936972, 28.00084236490239, 24.536704058241977, 25.478142140180623, 23.769089817522975, 23.636632846169448, 24.43084697534719, 22.53642327737772, 22.20977720091904, 22.03955648996786, 22.100537957126818, 22.902583198467383, 21.971717764043643, 22.957276402302774, 22.916387696102294]),
array([25.435646881641297, 27.076808322421606, 26.157481785645242, 27.65411692265945, 25.983077870819933, 24.883687922627935]),
array([28.651061334215953, 25.987028554076712, 23.43831469567826, 27.72073267087312, 25.7636104240937, 27.125364914568816, 26.952911129681194, 26.40496744835273, 24.510118665250754, 28.02164168456017, 27.968783868938814, 25.17258760288933, 25.2607051643409, 21.19165466035103, 22.427181289004963, 21.43264074216296, 23.00368215550316, 22.214270929435553, 22.635640704555282, 21.292439085368837, 22.883430542747035, 22.54074320019652, 20.771356873282276, 21.267024782884892, 21.41099002965024, 21.21783653444487, 22.130057127433695, 19.079146514898724, 20.039966483652844, 16.87724564164172, 16.21817484504216, 14.936352483641553, 14.823854959768113, 14.068223564710893, 15.06565374369653, 14.616474814565382, 15.59945844582194, 11.638863342934187, 12.504744180468432, 12.171087571837761, 12.213028840423117, 9.807491111019862, 9.207733945377262, 11.146053048892087, 11.55920483585987, 11.621455133853972, 7.868643578709785, 8.266281975682093, 11.284144346557365, 10.89083541491302, 4.972009979005661, 3.4567350659563094, 3.4654304411860863, 2.6738204330825415, 2.654601301676305, 3.2719544481162535, 0.32062646563012354, 0.00727457881429508, 0.0]),
array([26.480588779011107, 25.05913840993927, 25.841844480796098, 23.956762886433328, 23.476608566236177, 25.54760828116747, 23.075018191039163, 20.538254329679404, 22.48252338752442, 22.928354507073692, 20.693559610373278, 22.73171919500091, 22.95699459478439, 22.57545039603251, 20.948404324018682, 20.787575041964985, 21.769359053361462, 16.719723025015625, 19.999022556217017, 17.04272015830539, 19.65781470079409, 18.07059169538126, 17.986596219837054, 13.952618164069975, 14.009702410796423, 15.528026013133971, 14.573333011120944, 14.23717978308548, 14.348680708040199, 15.163023528891689, 14.79845064443124, 13.056328827623103, 12.140085637461013, 9.914602712709996, 11.356504833336325, 10.6164903110785]),
array([27.00467377126333, 22.40206320887879, 20.76210705313873, 21.16280014377942, 20.56592250642962, 20.605115991217495, 20.777792432800844, 19.13011317607651, 13.86461185750614, 15.849201648691588, 12.988471415558369, 10.909145956321684, 9.821929658440132, 8.112024569602957, 10.777725739052773, 10.106890996328197, 3.4886413027288237, 2.7192503992873034, 2.8608638274651477, 2.8537494104028354, 3.0987691364838623, 3.1061005955682894, 2.769191888762087, 1.9715354538890644, 1.279209534785248, 0.5935931655009675]),
array([24.212201970685964, 24.349930880823823, 25.726672241047698]),
array([22.860847697268518, 22.217637049710792, 17.934051472816797, 14.528571166887753, 14.099703375236118, 12.450323033903013, 13.745862420418803]),
array([23.95016674745388, 24.439318098188444, 23.5144406457105, 21.162892324258316, 22.1073799113966, 22.91294592560838, 22.27174420909975, 22.881957303087542, 21.17346558565333, 22.52514008067582, 19.095576036924655, 19.405678610013045, 18.77300850047274, 19.092926985517614]),
array([24.613680411329153, 24.019163903321676, 22.643979617480813, 21.190880615655118, 21.56498308345906, 20.72980533448308, 20.888645364746893, 22.27722728249881, 19.32180357377873]),
array([24.38194439394853, 23.632783540402162, 24.266522354998177, 24.379544151294706, 23.064940209134793, 23.757155868456504, 22.74042281839558, 20.657413234190493, 21.595138301117803, 21.985748725023804, 20.57606943707492, 21.363466301034542, 22.425360124909876, 22.180830543508996, 20.605222377419732, 21.21977866536713, 20.961988044357444, 20.48251766451102, 21.111717442397897, 21.002915476436808, 21.20544785824824, 21.565801414485936, 18.030345591229505, 19.906443162969154, 18.471884867770466, 17.61859711047415, 17.327481909501817, 16.82528257312106, 16.498041708849676, 17.919867654851245, 14.776826429745423, 15.243363739732974, 14.899558325592421, 15.688044910216655, 13.829662025428714, 15.189051893062299, 13.492892320349386, 13.749961202335927, 9.494230639815976, 8.449170243442738, 10.420908724232245, 9.820581380279236, 10.209349179929724, 7.856529324890834, 9.367211315998375, 10.52591052461274, 8.41980828464939, 5.401386274031632, 5.9746399861587545, 6.479361595915085, 3.913664033659699]),
array([22.054269599942206, 21.975061341061977, 22.948896405624353, 21.456425121151348, 22.3591284034809, 21.002999544507638, 22.47880393049387, 22.248667106654306, 22.59212475856553, 19.63214853646074, 18.198680288553646, 18.22279171966633, 19.85837352719021]),
array([23.06960356101775, 22.705646896893267, 22.331226739218454, 22.16123390698541, 22.842709344679776, 22.162454586992858, 22.17794815057576, 22.107943719887476, 22.39982329028959, 22.726054175009722]),
array([22.535591714285264, 21.401434893297917, 21.23353873428533, 22.463900575273016, 22.197551024534874, 21.941750728502278, 21.43649300679978, 21.172999896850847, 21.174015770755386, 22.349115068709576, 22.04602138617516, 21.802419501621593, 20.956911315127588, 21.174566481321943, 17.639358662314304, 16.882413310553257, 17.889631775587446, 20.023230560878115, 19.535957606034614, 17.623636978196537, 19.539605128696156, 16.864546948891366, 19.19604446144449, 17.00964928439609, 18.02938233804641, 17.044758365425377, 15.219734420107475, 14.784852365633245, 14.409266873788816, 14.708520550780644, 13.827865379393826, 15.197087971616181, 15.421731471370261, 14.19276398346696, 13.340237272082298, 12.507710746677997]),
array([21.10146524428695, 21.27811070128435, 21.89994461480547, 20.5933352613213, 21.014062789752145, 20.850742569698145, 20.957266151855116, 21.937371461913347, 21.426904153730455, 21.985538650084486, 20.555890583233715, 18.384801472740154, 18.11037867578917, 18.78061247484303, 19.709245031691704]),
array([20.918029194940814, 21.583365002213373]),
array([20.56766440409251, 21.6016177516479, 21.881699427763078, 21.35279754386056, 20.47485990251337, 18.16443631378864, 18.37558379446915, 17.246414306494852, 19.89694967512562, 18.662709626722414]),
array([18.180341355281982, 18.4293919743658, 19.283748189164427, 18.34295945066016, 14.236842758850791]),
array([17.456816483833354, 20.29743531042872, 16.76002702720488, 16.991104076877804, 15.969502382588669, 15.912764802994607, 15.943832583456228]),
array([19.63758441330487, 15.781774338623752, 15.303704085782961, 15.830815512740392, 15.877725313266675, 14.561633286723122, 13.768518145433912, 7.439411443510838, 8.139281344204072, 7.924368705406127, 7.2755767706909715, 10.483427209237622, 9.288464075975751, 9.584899888359223, 11.45016867599522, 5.9481612559884285]),
array([20.019400656678823, 19.864956534072007, 18.67777929702853, 18.65847404041397, 19.73346993281663]),
array([17.755856131702554, 17.14930611789372, 18.137725504002958, 19.099126992357068, 18.54052953734798, 19.04191501833852, 18.709774729056473, 17.393735802450273, 17.823465743840483, 14.591362311620003, 15.033535044158294, 14.195641263323212]),
array([17.010834247312214, 16.50508274739603, 14.017244455289568, 14.038120726939901, 14.32243328809842, 15.638150643066046, 12.729562225019997, 12.747853111351308, 13.552814677235085, 7.8491055538154235, 9.502857279488634, 6.599004815406987, 7.167860548778961, 5.552752455100837, 6.199235150286145, 3.390478377901778, 3.5286121406851807, 3.2403520118713893, 2.8552971860509073, 2.0601402032396825, 1.6634601933788253, 0.47691807603983327, 0.712616481103944, 0.6811198164899676, 0.3435283636866247, 0.0]),
array([16.09516272787952, 17.877672432435865, 15.744639948281836, 13.966817665614965, 12.040550247408184, 13.735020436831453, 11.81524141527765, 9.142981450142425, 10.734936090935221, 8.159597275582666, 9.373982876282435, 8.454629082612598]),
array([16.004180809518928, 14.290799391939036, 14.241785126076266, 15.659550120309321, 15.750795900621716, 15.59005703197764, 12.304381989487315, 12.483190615893468, 12.000967840955086, 13.0122680503996, 11.575901575219303, 10.320240415341189, 10.701633086699756, 9.261396107420198, 9.580941306558138, 10.954971805228832, 10.561274100643221]),
array([15.019945554266744, 15.871192561502676, 15.256091373890117, 14.798338036912238, 15.202841971745896, 13.74071972267498, 9.363759962819644, 8.757003757442662, 8.281421508148096, 10.113019330777945, 8.499945803129293, 9.764447517774443, 8.814102289409028, 6.944922074614069, 5.152113116687741, 5.055807711940296, 3.2351779234854368, 3.1886680409646164, 3.3769879180636257]),
array([16.846893901866565, 17.390160246640786]),
array([16.375069860835318, 16.3993467560787, 13.626379898453209, 12.674242549086735, 12.774246219403532, 7.946538566686888, 8.96313826528614, 8.67284977531136, 10.509303693821813, 11.338055096948949, 8.339341603656083, 8.023858229879409, 8.914051134401758, 6.43150637371027]),
array([16.81189078685009, 16.21530185442792, 15.694062554555119]),
array([16.798713238573686, 16.522381956787996, 16.418963465693196, 14.62320835730737, 15.363052811259223, 15.475036876870826, 14.630400556436307, 15.854114489730494, 15.960581870305955]),
array([14.32951418604025, 13.974404235466887, 14.48631612463096, 7.45951881506021, 10.570685322711453, 8.603649282920694, 9.575622250913941, 6.453061961333557, 4.077822262719933, 4.495080452514053, 3.42262197492728, 1.9244075837432209, 2.556535734672306, 0.540007869892935, 0.517344946021508, 0.12291011174203496, 0.0]),
array([16.484046805548786, 14.77904105981452, 14.615771697229986, 14.86080693987012, 15.448824419220504, 15.93226033115344, 15.008719322163056, 14.79879247827634]),
array([13.985046048566055, 14.828482228094721]),
array([16.176452726189623, 16.181496892502068, 16.452933420936823, 15.139417354412306, 14.288987763925917, 14.125460492627736, 13.993480861940682, 15.810624210994094, 15.365106592865553, 14.905148104877915, 14.013489784203895, 15.134870170392908, 15.041661611211511, 14.464962918945487, 14.926352254416985, 14.729982650418908, 12.221365241513837, 12.250046656059968, 12.470893779669808, 13.132298730352808, 10.977522730446614, 10.453939111630174, 11.22041749015612, 11.289853117377993, 11.35676500835433, 11.198650402802736, 10.083712728612243, 10.998011305684976, 9.791161600118258, 10.038347071003598]),
array([15.584751706118112, 14.616692199139493, 15.258101667974636, 13.913167550961274]),
array([15.518418302597018]),
array([14.964102512352452, 15.726897005589743, 10.27156453744663, 10.217119509018868, 11.2329897176845, 10.636389712737762, 8.785264890303809, 9.424920469953472, 7.9102596131076295, 7.706759119990192]),
array([14.499093156565275, 14.354380247753634, 15.382271778293354, 15.756637788165294, 15.447842888643693, 15.637014105545965, 15.36014001446877, 12.826974257557609, 12.429320478727858, 11.629559657256733, 10.42856468960306, 8.448236229539805, 9.789209676679109, 11.462906101955415, 10.673742720476422, 9.769524596230431, 11.444832808913665, 9.082780559729883, 10.864763389783013, 11.150649712533475, 11.427322199592595, 8.473541373862393]),
array([15.091898332061977, 14.593147780989442, 14.468535223081883, 14.661104123635985, 14.534542089985804]),
array([14.591677212005825, 9.104313134937069, 9.39275908448737, 8.298699671912384, 9.146784164119792, 6.104762882995173]),
array([14.514598239013456, 13.519154922466305]),
array([14.42775476826332, 14.684901504199146, 13.24849152168229]),
array([14.227485335683923]),
array([13.872240263868914, 11.070498427402404, 11.282638708760265]),
array([13.844482855441935, 12.930714691694927, 12.340619663568281, 11.639250089453082, 10.901675580819528, 8.058245465336473, 8.140017442259373, 10.80979681657961, 10.321763849405869, 10.295254509304465, 9.482620318862034, 7.882799561056068, 9.64701525881678, 7.291015550262823, 8.976739969786662, 7.588388062930459, 9.44502319953987, 8.229487004978985, 8.604606671633071, 5.3482229762971265, 5.744353964971481, 3.458964747237452, 3.4589938775469298, 3.2176047454825962, 2.701028235118485, 3.33765603370169, 2.0807514569870467, 2.5615857224622447, 1.9568882621588706, 0.8308996586455143, 0.6838586645058158, 0.5198606895566802, 0.03231482705935797, 0.0]),
array([8.455335030215302, 8.308206958852624, 6.135674700699646, 3.3402840271695227, 2.8000865456686705, 0.0]),
array([12.303530200559878, 11.724510841679058, 13.083454849731249]),
array([13.113282085494276, 12.7285305750724, 12.530327352711845, 10.646402980061753, 7.917927213334417, 9.675524233479088, 11.49998507494955, 6.366892634771228, 2.8064300337605417, 2.6220440644806713, 1.8545135442383174, 2.209362349428483, 2.425145634153084, 0.6451554189910577, 0.45309294610610557, 0.0]),
array([12.973968548960839, 8.320107858731213, 9.559767231035037, 10.832332145396899, 10.317882103627493, 9.36561326590667, 8.815856989488392, 8.53093969891154, 9.506844786962182, 8.644674784246714, 11.368319879164707, 11.035318590576622, 7.122156151075748, 4.516942298087931, 5.132612042017497, 3.0581783673466636, 2.6343136842177697, 2.2783536690355293, 0.8773466377051183, 0.23363890233704088, 0.5674431517978511, 0.711028090156532, 0.6529170662222618, 0.03056805416237418, 0.0]),
array([12.452240400724596, 11.244152090851843]),
array([12.908629360498425, 11.680122626428894, 10.43786758833172, 10.718431069204213]),
array([11.73855910394672, 9.331288067032194, 7.2679111196486845, 11.468176024249074, 10.055655821872502, 10.162539196484042, 10.513461768372025, 8.965122129256983, 7.247401615900257, 5.769016916008881, 5.9916262152644615, 6.784745672327555, 5.586374278480817, 3.7917910564842914, 2.7767084978988397, 3.5618044266826354, 3.30007112141855, 2.8169691102576033, 2.696955364451369, 3.0265171085868, 2.8222726669400546, 2.031516288780211, 2.573978831994153, 1.091593600546038, 0.22907954677018716, 0.3632946368022837, 0.5007066005087499, 0.3373019140539062, 0.4063646454325605, 0.05423941390103207, 0.0]),
array([12.406361169851099]),
array([8.676530917492407, 9.85902295366417, 8.354760576277151, 10.976561862642075, 6.393581234291761, 0.0]),
array([7.346762582389761, 11.497239825904357, 7.347449553995139, 10.632657453407093, 2.6962668913108625, 0.3447139304519891, 0.0]),
array([12.44545031514091, 11.159791624343892, 10.55378215253057]),
array([12.136718714920725, 11.170236917761112]),
array([12.042716151067163, 10.512615251380154, 9.986039748902451, 6.258491131458142, 5.242261660776977, 4.47505170587708, 2.8901940536102875, 3.3039550117057233, 1.9010928831966813, 1.5740460551659725, 0.5551023603293168, 0.0]),
array([12.100623005812608, 10.750841489375272, 8.944099621942149, 9.540625768618376, 8.02400045233508, 10.646322227703047, 9.417564316087711, 8.463100791499524, 7.831523447278052, 10.132414342541345, 7.926489804646199, 10.996752501202836, 9.66672192784066, 10.223272432961322, 7.676532911174508, 7.969300893456257, 6.436745808057767, 5.982081264043224, 6.018053455397389, 4.060256528983627, 2.89775471710678, 3.4574528547803722, 3.563324640936266, 2.683714880510067, 2.2469905452373347, 2.033009089060323, 2.16973468616132, 1.980312172063706, 2.172479350596776, 0.8873666148240932, 1.1919334391903083, 0.7430309905339949, 0.4406461818359801, 0.19025118467511604, 0.24226268615521218, 0.36535495055710526, 0.4878549034477743, 0.6429470327677407, 0.7154465917177075, 0.6810475878150162, 0.5857467552965167, 0.0]),
array([8.813888131633497, 11.450044373576727, 7.167596285998444, 5.514663300353966, 5.321643532099699, 3.3249750121066115, 2.898087514494769, 2.8130493593020933, 3.099953923226086, 3.197837362144356, 2.3195746938948476, 0.49310208541205713, 0.36323157330396827, 0.0]),
array([11.500102108840714, 11.179722546558708, 11.120710115796914]),
array([10.541034944964276, 8.931945441069974]),
array([8.643410613550522, 10.031998521178991, 9.636836022401818, 10.611945512857291, 8.067708491962211, 8.896924229620055, 8.596587295844394]),
array([9.437461456168789, 10.254746004361381, 7.627540144576335]),
array([8.732220332171726, 9.636518641953796, 7.418706759962904, 8.967615397195871, 9.036327648439743, 8.220996791879827, 8.12877171599464, 9.918380407820486, 8.384617849252528, 8.26204461444968, 7.688400373488159, 5.660488043480866, 3.332501234108387, 3.47513904096148, 3.5386952819665694, 3.0802942397962942, 1.9293353105986073, 2.4213524838436253]),
array([9.883204659670126, 9.755487672577049, 9.379272002211525, 9.673122704291607, 8.963455385325737, 9.789434626024553]),
array([9.129002267525276]),
array([8.550384420290632, 9.735818937552057, 8.157186296979155, 8.983176814776403, 7.215823139558608, 5.0323361042505095, 3.0126930895621333, 3.1195241531104454, 3.391751914029861, 3.3468615091464238, 3.2398031358588515, 3.0015717254201695, 1.8262662984373685, 1.9617648022442498, 1.1470124649567524, 0.7690652866497839, 0.5100582855509261, 0.335875767920045, 0.25153998938680666, 0.11047822158949415, 0.0]),
array([7.952869052942857, 8.327287295045029, 8.942899727421848, 8.463380931944336, 7.978554864920108, 7.572603064645721, 8.688340456331428, 8.419994357798881, 9.12537754495017, 8.799450265799955, 5.893046961000624, 7.214635548023428, 4.886776197269272, 3.75442456635826, 5.246753740380807, 3.2529679294943272, 2.7140700784540583, 3.0518380769634126, 2.395793481919902, 2.3334738276353266, 2.2149948118402354, 2.5607148332681806, 2.542137363984645, 1.6767421970947476, 0.20462496837331712, 0.23280830000554187, 0.633452067130076, 0.2032830522607455, 0.6442235904675894, 0.0738534458276543, 0.09056536315120198, 0.0]),
array([8.552460418100987, 8.385359865582458, 9.377526904945695, 7.474802341025484, 8.152611388622228, 7.745346700468902, 7.728486957307495, 9.080941253678796, 9.096042215684179, 8.053476009627662, 8.370976386458837, 8.918110441579493, 5.912461679833686, 6.598972659635835, 6.991886695717636, 5.818971073962784, 5.437951750672003, 4.300467547755904, 3.9633883126082705, 3.485889772046268, 3.160218880060809, 3.053037144987562, 3.418319346816125, 2.937262655700779, 2.2927706137336514, 1.9073725387003, 2.5110899393617925, 1.0968850407726207, 0.6697818590302156, 0.31067795227692124, 0.5459043336766913, 0.45271169397134037, 0.3439605122016103, 0.6157483676105313, 0.6451848162250631, 0.5356219806087411, 0.5646426560841231, 0.12002543280068596, 0.0]),
array([7.259758813736633, 9.557841345879801, 8.8585944157836, 7.5821128921711285]),
array([8.3539190674044, 8.349365881180711, 9.281985935870605, 7.918035559568233, 9.254389143868249, 8.391017465010835, 8.903978745215761, 7.9058885404630965, 7.424824903028719, 7.757280904790585, 8.43873432933352, 6.922417935207311, 5.601952280591674, 6.148581611461018, 6.382468546740468]),
array([9.109590515574263, 8.911490581006834, 8.820987300355288, 6.271400756338477, 3.767935653494873, 3.44609025902958, 3.270819489358099, 3.390485441981842, 3.3754854540388934, 2.31712486374609, 0.709775486060623, 0.0]),
array([8.833393154119834, 7.487157353496458, 9.158682802267489, 6.636845064314415, 3.752752073908023, 3.3938838048613493, 2.657169687585135, 2.909015724786088, 3.5964846356145106, 2.848160861677731, 3.299995735644238, 1.2660477935357144, 1.0710576407012933, 0.3922843868364524, 0.4316680403778987, 0.45182905037064514, 0.5163859291966884, 0.25222662918435557, 0.19006263523794598, 0.2066501064822669, 0.0]),
array([8.734267060110032, 8.55430492415564, 8.381680623098461, 7.510000363982445, 8.083737616128634, 6.636877982794666, 5.870435856326414, 6.250173504643223, 4.046338036968161, 5.148834970151885, 4.96572234398919, 3.451531717518249]),
array([7.472315800801821, 5.563306461276445, 5.556832913386899, 3.7572525147780658, 4.105892105662574, 2.8721145779567747, 3.3067502792873675, 2.8334097522309616, 2.1804757448953924, 2.576777455874613, 0.5395986853392779, 0.06822139956462547, 0.0]),
array([7.6316157403785585, 7.340234654962056, 7.747270759918329, 7.498412512222436, 7.892665702302493, 8.362148570512174, 5.7684249581724725, 6.5318939258081405, 4.323671930042814, 4.592066666970126, 4.820535181929274, 4.8904394904970205, 4.457654262125969, 3.429693036513175, 3.186277351195181, 3.037328909031223, 2.5982477634491925, 3.4296955578227895, 0.6287231623583042, 0.5746599927797461, 0.0023295768887159757, 0.0]),
array([8.007764180011952, 7.608220040606653, 7.771372135786325, 7.419432047021483, 8.383909701198377, 5.910523537979268, 6.724672591097737, 6.331168755140856, 5.035261346278686, 5.142467786657032]),
array([7.945973133405466, 7.473163355984698, 7.512394141887549, 8.539545913076976, 6.869781803750773, 5.722320954274494, 4.2901454158298575, 4.163205344954191, 2.9365572450391157, 3.250038599098277, 3.2305644157421014, 2.9687349788165935, 3.1083515298868845, 3.0461128448669363, 0.9684360729001047, 1.7589634996659713, 0.6162276227665436, 0.7352001814417637, 0.6659944770505163, 0.44075018374169467, 0.48739153064989704, 0.4685243046702753, 0.7198521382297014, 0.6756530103046958, 0.0]),
array([2.7015346912160094, 3.4772665494367185, 3.22401667206601, 2.4679105088855304, 0.9475826997938758, 0.0]),
array([8.14535414768509, 7.472867344295355, 8.11441361731359, 8.135797918482698, 7.512115590130291, 6.398881340322533, 5.629041845748828, 7.186446143228569, 6.669984908533536, 6.826219412645936, 6.448916122941159, 3.9871094047805506, 4.950663382172533, 3.7089048667237803, 3.5908333404684227, 2.66310445223638, 3.519699590010171, 3.532481940554262, 2.6879400791027117, 3.387302005935657, 2.846112788894988, 2.8659956260137753, 3.333388832324834, 2.819930398313371, 3.034144330396101, 2.8877753247175195, 2.3546325073381484, 2.418588003857356, 2.3626122739391717, 1.1769921373768355, 1.5788582767441401, 0.5351004890390672, 0.699560850628277, 0.7400285971379896, 0.18344632044520826, 0.15093341739361466, 0.06775386402678217, 0.0]),
array([8.242710827702156, 7.4944853532608136, 7.973568982773162, 7.457787757987966, 5.403015948039246, 5.897413526398545, 5.539190426310661, 4.631096530400284, 4.059499472809834, 3.129738989369373, 2.7389520789719173, 2.9463682261886452, 3.487235861257844, 2.84841619616867, 3.459745143760122, 3.330493919275719, 3.0919205133266683, 3.1029245468571625, 3.06037885214132, 2.751949742663042, 3.217124978323212, 2.5135932761199093, 2.412751598710594, 1.1445187438042435, 0.16171356552316962, 0.3359772110502996, 0.30497713465705484, 0.503806533721542, 0.35684911923354623, 0.6389001016258826, 0.32685410610360277, 0.6094821730269844, 0.23044032472043996, 0.13286502700181646, 0.0]),
array([7.787842956793048, 6.2703600069445145, 4.307493846952874, 3.88031029800658, 3.1532315391822108, 2.969860973978908, 2.9123501740625737]),
array([6.401367464036066, 5.403879531794888, 6.393203735506242]),
array([7.238275208995597, 4.078393749143281, 3.141581547877495, 3.476744101102108, 2.7503070830242824, 3.4127613461497663, 3.3780339249081086, 1.8115788497489294, 2.354088121758582, 2.5728195708488233, 0.7342305842790733, 0.49561863261898786, 0.10463900966562352, 0.0]),
array([5.7037830508123415, 5.53524451498631, 4.216561074781136, 4.270133210096963, 2.8181643296983485, 2.9016292598213482, 2.8118117795322846, 3.2904533712847552, 3.2131848136658494, 2.6751632085205, 3.489215317711192, 2.768375715118956, 2.789118704938307, 2.165471158005449, 1.822805449698114, 2.034051103487279, 2.0530779119044786, 1.5420143861336553]),
array([6.072157136417568, 6.091653415406821, 7.061276916933861]),
array([6.4260934470990065, 1.8881901133170094, 2.042740674581437, 2.1293348257882516, 0.5083722028558832, 0.5268845937082386, 0.7403239628343723, 0.0]),
array([6.4209370329935584, 4.871825663263787, 5.257067470131465, 3.2794316192236943, 3.5185723895385896, 3.1587550672119504, 2.5601533051365597, 1.5570639816722234]),
array([5.565801812123152]),
array([5.791944354962318, 6.826054260236332, 6.865555455577772, 5.546638076966344, 5.87054623644293, 4.624749800346144, 4.914737737526128, 4.743992971810179, 3.426930038509452, 2.7532900634870865, 2.6065108699815163, 3.4020704764294174, 2.8992054925580466, 3.001503280613586, 2.9833069869059416, 2.970928723585688, 3.214117448218503, 3.152840194182426, 3.3013901842579894, 3.1130697230611886, 3.548779369418335, 2.8189813300604802, 3.4106056256825568, 3.5187375285813207, 2.0089350081057376, 2.4023429955870137, 2.2953884504613815, 2.331397386857124, 1.891234010834133, 2.0687746102211166, 2.506271194726447]),
array([5.786816983090108, 4.251490521309055, 5.246788838420598, 4.5693662157481425, 3.3591783833782065, 3.063552397953429, 1.9425629393996053, 2.257072105952661, 0.42481322819755984, 0.15454100020303607, 0.17233819998265787, 0.2318852448004396, 0.0]),
array([5.024788035758052, 4.0401033647000695, 2.691320834448563, 3.329202805671543, 2.7882255167057375, 3.202829019458051, 3.3759337411436907, 2.750280340759904, 2.648775944035047, 1.7336889049305657, 0.29852409540018904, 0.17128215116480872, 0.15321611156622073, 0.6968853917091896, 0.7066740278365745, 0.0]),
array([6.32359459053346, 3.020848855528789]),
array([5.5432640026152145, 4.42319746351416, 3.0130536474848633, 2.7246860344329082, 2.1868642475382085, 2.49064167305557, 2.351965191542468, 1.4071208390412937, 0.28017152632172304, 0.6140958949922527, 0.25975393030856986, 0.0]),
array([2.9321749394240273, 3.1934710507317843]),
array([4.622458186443266, 3.589475573937319]),
array([3.8687159400619446, 3.3100012494504334, 0.08122773174525588, 0.0]),
array([5.429250465641609, 3.567832339960495, 2.986818377266191, 3.339882852930122, 2.8662147281694903, 3.0699735547296827, 3.018193982925039, 3.351611201490216, 2.6642250996194545, 2.936569242310632, 3.5141489512691373, 2.033316716498337, 2.556289730572291, 0.5100668530340264, 0.34355032467585506, 0.72686326717013, 0.2762399171847889, 0.226631702464682, 0.39549434475247663, 0.08080067986162184, 0.0]),
array([5.358374166755931, 3.066637910195408, 2.969862578935744, 2.5680058211373322, 2.574916612584517, 0.0]),
array([0.5676314939550003, 0.0]),
array([5.147603557613463, 2.915040518178358, 2.5805845017282376]),
array([4.465493782528328, 4.239709040315637, 2.94698306572737, 2.8908037258005668, 2.5730955229628525, 2.376814911926253, 0.6955264174850389, 0.13370295685005296, 0.5417920519481261, 0.510990263120578, 0.5138411020760394, 0.3156283180244961, 0.0]),
array([3.7449703179119087, 4.220470073593574, 3.0427454871869495, 3.111788606646752, 3.50911453608104, 2.8650793295830477, 2.188522232346804, 2.1318698911355503]),
array([2.8434025415736452, 3.0586151647432174, 3.156996179431172, 2.7244703758632225, 3.4387272900614096, 3.354118869461126, 3.5802177248371954, 2.4300101777243435, 1.8558068147641071, 0.7222119532573684, 0.21194709081497176, 0.5554696472448976, 0.32172023527605165, 0.024044885328143098, 0.06612046877247073, 0.0]),
array([2.9605193265735887, 3.501622362231634, 0.0]),
array([3.433762827069734, 3.2900180604300857, 3.339370119432425, 2.8255373898614997, 3.4675263494317607, 2.7791005808917832, 3.475189423038633, 2.6621386341827806, 3.3778277476574443, 2.2405945264594207, 2.265510437952653, 2.113614510699306]),
array([4.697057311757769, 4.22685139319319, 2.8191346648757762, 2.8259927736189328, 3.4846602803370708, 3.1912922862515725, 3.314792878396863, 3.368395571377127, 2.322837932634494, 2.037822139695516, 0.9871167384687625, 0.3008534754380301, 0.46689124451215314, 0.7545448203156204, 0.5419840406249359, 0.033430784233616465, 0.030823334250118142, 0.0]),
array([2.6910053633496624, 2.9885912466103557, 2.622671654014601, 1.3414438516065346, 0.1770625074822193, 0.6261112705556346, 0.7211123461563838]),
array([0.21401567288904888, 0.6955864002251138, 0.0]),
array([3.727901689456008, 3.707813785246665, 3.9106732635415966, 2.9284013376598894, 2.9595807390928752, 3.3406890477425777, 3.174648181096126, 3.531968893649065, 2.88535828253993, 3.270728260300302, 2.1727527697892417, 2.393873758365924, 2.5232278236071237, 1.5140293242255225, 0.7365992728490419, 0.20872528043988225, 0.3659292677444036, 0.0]),
array([4.284710449102899, 2.954361964388063, 3.026256460073134, 3.1887170054081886, 0.4307118596332767, 0.5306763182778487, 0.17478836401025277, 0.3989375813759768, 0.39528824440054583, 0.0]),
array([2.925600522572861, 2.2353109257781805]),
array([2.7862906200799906, 3.5139803274481887, 3.3440026804992877, 3.291790557999301, 2.6450825645497806, 3.4271005615594934, 3.2417850859367974, 2.872048462906932, 2.348434564929402, 2.5251922594528016, 1.3048438596732062, 0.6426641094053236, 0.1922780175600627, 0.1458144135595808, 0.14241031919551084, 0.33353319271939486, 0.2045225144684687, 0.21163056935624303, 0.3335956506178704, 0.0]),
array([2.6720642004210724, 2.489353388062804, 1.4125401304933802, 0.38193571558189576]),
array([2.9803732983473337, 3.1259073355319, 1.8312435929418038, 1.005150388302626, 0.1791903807251941, 0.6609232821448854, 0.0]),
array([3.8553440271832593, 3.6178874862080876, 3.1799133071274217, 3.3502488976062827, 3.3653970462171987, 3.319618213538911]),
array([3.314147472490637, 0.5281091683604464, 0.36212789377574445, 0.6786289596973246, 0.0]),
array([2.791218654329402, 2.775150939254418, 2.2466378858124827, 0.4363872410372464, 0.0]),
array([3.282369733368291, 1.97480806962988, 1.888450428789774, 0.39682076931291826, 0.2531637658320858, 0.7591908471937213, 0.2812751466155141, 0.7602062783318153, 0.0]),
array([1.8100907998879394, 0.7678758183933764]),
array([3.4147950896762875, 2.8753558504451266, 2.921595947668537, 3.224471021141792, 2.1190611136796984, 2.024391439284104, 1.5433756428511065, 0.9787535935973448, 0.7332703918117294, 0.7806293631477131, 0.2442295055163336, 0.18903587143474287, 0.14847580179725106, 0.4431185117899236, 0.0]),
array([2.753698107511564, 2.401552142934223, 1.8985168414156846, 0.5315263925415284, 0.2889216948091199, 0.0]),
array([1.788849933804345, 0.7341337414544912, 0.6255329751434192, 0.0]),
array([0.4061248420609535, 0.0]),
array([2.8944776328684303, 2.154356357659829, 0.6437976827399197, 0.0]),
array([0.6471944494655381, 0.0]),
array([1.9578274031606808, 2.4007863264313887, 1.2390381246657973, 0.5940981303246287, 0.6096642368192075, 0.3858370174531545, 0.5087332868132535, 0.023810278276240646, 0.0]),
array([2.381525441190978, 1.230923285673136, 0.32729230643154694, 0.43889640405052943, 0.5756158955376343, 0.26304035895667166, 0.6597462467694053, 0.4863452481146595, 0.20962122609934852, 0.11248814565781667, 0.03656542788560989, 0.08651717199616801, 0.0]),
array([2.611764568201766, 2.6299099035303883, 2.6520805684346738, 2.4541800542549863, 2.1953784148874744, 2.3050037795602334, 0.18656369050713428, 0.21441970987132752, 0.22006507377427464, 0.29877019620476025, 0.5087933923055712, 0.061906582970480156, 0.0]),
array([2.3443698487533533, 2.130763770375838, 1.9704450587062359, 1.3965592366749924, 0.14923072355369438, 0.20743988432428928, 0.0]),
array([0.60923357632567, 0.0]),
array([0.7548978233889357, 0.21214126517808995, 0.5088923606333561, 0.23027803949433934, 0.0]),
array([1.917446180220602, 2.008731466544986]),
array([2.0028299686169984, 1.8369104184645855, 1.5876625179463084, 0.22758175964554195, 0.4718133291268896, 0.7121825637379757, 0.414308936113081, 0.0]),
array([1.0843465090434243, 0.4354050965585816, 0.0]),
array([0.42446483855462974, 0.4487786448820832, 0.0]),
array([0.9843293167046504, 0.32272370230186026, 0.18516740812085786, 0.15557254659055553, 0.47062489302198185, 0.14366849825455552, 0.6305470906773041, 0.24947126391836627, 0.15420410867186074, 0.5843424471997749, 0.24772419826986625, 0.0]),
array([0.7581227528348079, 0.5426653906167683, 0.727666035459172, 0.0]),
array([1.3623113933124353, 0.4602920370704281, 0.2628063294144113, 0.31292781947583576, 0.6089858740101592, 0.7341945143008415, 0.29369928984646765, 0.1733261132853514, 0.5235408071012551, 0.244913604427337, 0.0]),
array([0.490074812921379, 0.0]),
array([0.20068961356177006, 0.15714811867562972, 0.31298218860817684, 0.5206108598890293, 0.5604794835223972, 0.6327593084358476, 0.14381062277528223, 0.5659253940970757, 0.0]),
array([1.134347961478018, 0.5488740185412562, 0.7649402493213711, 0.7058686420808429, 0.486691322219338, 0.2515317324372811, 0.10000930174740252, 0.07744679380357258, 0.0]),
array([0.5515135671467394, 0.5016812538903674, 0.6166884538590757, 0.4267085615987894, 0.2757521218120793, 0.4849559831182361, 0.5607758929234796, 0.7208591405692643, 0.08945693750812823, 0.0624264237472832, 0.0]),
array([0.14414780696414764, 0.7150923439862631, 0.46716770414777703, 0.0]),
array([1.5215162570860823, 0.7771882739229954, 0.25679773964050145, 0.0]),
array([0.7256773621805784, 0.2493398380774834, 0.6661411804674752, 0.7590476789134666, 0.0]),
array([1.3030649413827504, 0.5171687222754505, 0.5168771777796286, 0.6694723489544985, 0.24365398130631744, 0.43922396734599706, 0.0]),
array([0.16264642685922026, 0.22589840063158195, 0.0]),
array([0.5091246583636317, 0.0]),
array([0.16658974371567914, 0.21642082145559083, 0.0]),
array([0.4800765341960582, 0.2042289829649615, 0.3280878994205863, 0.2101391085432811, 0.0]),
array([1.0654598330293905, 0.5435876397825383, 0.3052712510094974, 0.7512701470054848, 0.3228466730242765, 0.27574225720577517, 0.6803261530756509, 0.19926186580119953, 0.19633588466484297, 0.0]),
array([0.9236869372016598, 0.6932815821849747, 0.0]),
array([0.18287782262896446, 0.36843939356538613, 0.1530864699651302, 0.0]),
array([0.3607153270109861, 0.1509075505111348, 0.5222824507340148, 0.0]),
array([0.6123030290990732, 0.42561271320593885, 0.43916016167407623, 0.1184078634200406, 0.0]),
array([0.35025781692681635, 0.2508167259339551, 0.17122393055427687, 0.6072648944955346, 0.5256504826953472, 0.2794388858171315, 0.05385817888806041, 0.0]),
array([0.2672540729481234, 0.5722197797273385, 0.37429341418006634, 0.08606191450164735]),
array([0.6789706490744507, 0.021822924590895243, 0.0]),
array([0.4391097291143371, 0.5410792630593728, 0.7573732608491478, 0.7492472866165082, 0.5898764420447687, 0.35344940165779853, 0.0]),
array([0.36369881223602074, 0.2187388482208017, 0.536937203649255, 0.662107028583299, 0.0]),
array([0.23604947091758538, 0.5712048864695454, 0.5690335406234253, 0.17869081587107472, 0.507735019470875, 0.07364739290348409, 0.06597483088258363, 0.0]),
array([0.2901232493448539, 0.0]),
array([0.4520085080271966, 0.5724128808840377, 0.5748277188422445, 0.4794939561955023, 0.35578381748217464, 0.0]),
array([0.23692576269107618, 0.0]),
array([0.3612423271019201, 0.21856980449751207, 0.009373257657128348, 0.0]),
array([0.3066715933595823, 0.1429529153167357, 0.27211426078838546, 0.35520927586523093, 0.07144444802594921, 0.0]),
array([0.3458123213082521, 0.36772181743871846, 0.3200313180943824, 0.0])
]
d = [data_1]
names = ["62"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T5', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T28', 'T29', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T43', 'T45', 'T46', 'T47', 'T48', 'T49', 'T50', 'T51', 'T52', 'T53', 'T55', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T63', 'T64', 'T65', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T130', 'T131', 'T132', 'T133', 'T135', 'T136', 'T137', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T169', 'T170', 'T171', 'T173', 'T174', 'T175', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T188', 'T190', 'T192', 'T193', 'T195', 'T196', 'T197']
def get_taxa_names(): return taxa_names