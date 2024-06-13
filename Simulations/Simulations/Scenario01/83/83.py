#!/usr/bin/env python
from numpy import *
data_1 = [
array([32.51840225886285, 31.44278495288031, 29.322759508837805, 24.38959901750547, 26.899381645776568, 27.17915646643371, 25.625919439116394, 27.472171389860424, 25.610294474328963, 25.23719845255636, 23.69598421929353, 22.219892299906753, 21.205669001923678, 23.026617440219532, 20.764889556885674, 19.38372964369195, 20.347170428034282, 18.72117436743893, 16.45838027515712, 18.998107683825406, 19.431891542103294, 17.27565647326449, 18.101216690185712, 17.395471633436113, 19.25523624673754, 15.941095077006768, 15.687406680003406, 15.885465827116901, 12.499944932335792, 12.47715146560208, 11.006645599909673]),
array([24.906499491730667, 25.780253610560152, 24.7241336346641, 27.30044367774762, 13.897379362283182, 15.06705236183271, 14.253933490288171, 8.191057561893071, 7.519877786501677, 4.415317231241427, 2.080186280877733, 0.0]),
array([27.489269476289145, 24.70885443932242, 25.546500330226735, 20.587438137028087, 22.541554610624743, 20.9496906268268, 19.79215631837307, 18.643965344390885, 16.954528590216455, 14.132715280484437, 13.58285921459243, 10.346124071595556, 10.218315440783488, 7.437786060801946, 3.5764728236226526]),
array([26.300909328294924, 26.974095360426592, 26.830344714939503, 26.06434841404013, 27.78960455430848, 27.8344244302188, 26.230538792767987, 27.818363019592358, 27.164255794244223]),
array([23.75746519003031, 24.36372309135956, 25.801726056963375, 25.155291359861923, 24.265344067269822, 23.335813703531716, 25.23536017700415, 25.09749940498114, 23.412406382237435, 25.64006271614302, 25.130865650411813, 23.428434961046868, 25.360526216715794]),
array([24.297233863750478, 24.057326081674525, 21.59178881083748, 21.45107977658139]),
array([23.751915306287913]),
array([22.62315827218542, 22.451873658668873, 21.031464316235564, 16.632669875127863, 20.1733440218549, 19.41616913996064, 19.212614148033317, 20.132407865293516, 18.30012220079045, 16.73481772647474, 17.980433858498323, 17.975635880219265, 15.360581649161132]),
array([23.112168843179113, 23.661423124198514, 23.74533011977175, 21.679781706399712, 20.44086462628242, 20.98921396243283, 18.459396160663903, 20.21082680174504, 18.87562658944135, 19.368746798169873, 17.74950591819202, 18.14350129027524]),
array([23.14505576981609, 23.731847166985933, 23.178810318325855, 22.95288459061112]),
array([21.47457417337549, 22.64061497132122, 17.26029054295049, 19.80150893107988, 19.00215972717054, 19.974101044269954, 17.940285696338492, 16.3483551384827, 18.611798470449088, 18.92418106325068, 20.2065252230381, 19.982171201237463, 15.541953134783766, 13.873467882527411, 14.626126161822823, 15.195921726954568, 14.536173976734464, 15.64706936713293, 14.684042765559624, 12.959145442447753, 10.372313102376946, 10.625920099017975]),
array([22.82899626897691, 22.61644458004879, 22.673052569659532, 22.203515813170277, 22.35203552049435, 22.535480699071552, 22.56338414490909]),
array([22.021650925684046, 21.362109666718897, 22.104368718806477, 21.170505784119232, 20.932896016961895]),
array([20.699912749560465, 20.99736555377547, 20.98964707938637, 21.146706239937313, 21.255663388992225, 20.652188209655442, 21.52646818826105, 17.679249576611323, 19.821867017679057, 18.109667684071482, 19.937412922678952, 19.948200474940005, 19.05817160314868, 18.439767133171976, 16.155451262222783, 19.23720203784599, 16.476929702368576, 16.707441417355046, 16.76600282101037, 18.474866655153615, 16.352720256738905, 17.751624712474786, 18.73947816197006, 17.372311565924427, 17.20364550797946, 20.346632471165858, 18.48676399810559, 17.801903270960448]),
array([16.093476039471412, 16.761097177251198, 16.382254234420202, 18.94913777058933, 17.452995148760056, 16.857379340539655, 19.223420380993545, 17.714496477331835, 14.008521942316763, 14.256987880589273, 15.874385806192306, 14.995579575908016, 14.995607963271151, 14.094566697744476, 13.933383140386484, 14.260897245059082, 13.977103149567842, 15.637283838078332, 15.07239376822079, 15.224138394339963, 14.18673170497119, 14.620151636674112]),
array([17.143504192186075, 19.292539824567612, 17.63886451559742, 19.207651975992253, 16.026704947277956, 17.51942946600675, 16.726294281086737, 18.40699883798574, 17.135936795621365, 17.03516070182302, 14.257924382683303, 13.964745675400714, 14.817634041201293, 15.889345266554434, 15.638738617116225, 14.61477992918194, 11.671040665216488, 12.384994025696662, 5.619551438980739, 6.348358615352784, 5.872759327804521]),
array([17.276838098391984, 18.928810109274995, 17.7344745569031, 18.484889438510034, 18.407246807370136, 17.46944471242613, 18.23464165525304, 17.9041725688172]),
array([18.628807178533954, 18.221257067371212, 15.12328555523464, 14.300151930982164, 15.59581317064879, 12.793415139769502, 6.351343011290203, 6.061855072273212, 3.564597055420558, 2.478363513683776, 0.0]),
array([18.17426390579892]),
array([16.603417308222102, 17.413875594690097, 16.60897384513333, 16.301267564719886, 16.732488047963386, 16.240087826358412, 17.64249496948297, 17.39646171022464, 18.107005787471877, 17.863786793341685, 17.090636588044585, 16.43265239980278, 13.843329813487427, 14.942867006103468, 15.320958553412503, 14.165733905289255, 14.907466906061543, 14.451841585706896, 15.304067058328478, 13.986539734832778, 14.339286393855767, 14.9122484526075, 14.216497584221965, 13.059103118499399, 13.696296868338365, 11.953300565442483, 12.848554326004827, 13.666848075387902, 8.352467424415295, 11.385326857321184, 8.763857835697799, 9.156201437910013, 8.123551537310044, 8.822419105859462, 8.023661029461621]),
array([16.092203979605767, 17.19159070386705, 17.37399503661452, 17.200363258999186, 16.037414305804646, 17.8361000578158, 16.384614463666917, 18.19194975285461, 17.18875623313444, 18.290312696928357, 18.203675435219214, 16.516361851902637, 17.91132539300502, 18.33578236164818, 14.528876545069538, 14.644166700583748, 15.647809824947936, 14.289032948227561, 15.126450100081003, 15.110911874616116, 15.39163593293997, 15.056379996129854, 12.750200206569497, 13.563468038820403]),
array([17.148569460187176, 16.47039107082597, 15.194806867114426, 15.562771734963812, 15.068958860980892, 14.454827324039377, 14.158116816458026, 15.598164398265892, 15.169648778682282, 14.419208451613873, 13.5539877566046]),
array([17.497231134957833, 18.00162696349219, 17.512521541610298, 15.98973953248781, 17.24917017304751, 16.641085278966116, 16.855842048655887, 17.8915414536459, 18.159171800124117, 16.03575681771974, 15.610094364264548, 15.257999222034128, 15.629763029815194]),
array([17.915112758717264, 16.750675502639687, 16.175286408358513, 16.002734279534653, 17.0943781071468, 16.091310828398953, 15.647688908195686, 15.697261492036253, 15.322958509236853, 15.310100666592072, 12.998139858452127, 13.523214866517284, 8.78625294734262, 6.947612030707728, 5.964684359358181, 5.578633364122469, 6.817260977046301, 4.867363811950729, 2.5865965244125384, 2.9869724838499208, 0.0]),
array([16.070270168281763, 16.182571425389074, 16.86772776129668, 17.46596265654921, 16.409560528347615, 14.35443560586334, 14.126226904377639, 15.444762166698528, 15.430989487580455, 15.823904290080998, 12.545702222939523]),
array([16.47453845968418, 16.503512847743426, 16.965445895519494, 17.468174333863193, 17.067168611679183, 17.1928038623949, 17.79065768514608, 16.590473441830632, 16.624645317635704, 16.61989603300008, 17.71026991993501, 17.736360874311924, 16.99980240393343, 16.797199436002945, 17.63039617856864, 16.45278923863054, 17.553694610900028]),
array([16.341100488531904, 17.658486221393616, 17.674923166254608, 17.270155691588116, 14.502000758620976, 15.337992947615348, 14.814726129197716, 15.104237233015832, 15.20719788749818, 15.909710222308755, 15.162061158186907, 11.841636195680564, 11.719014180416057, 11.159777409480531, 8.820598320031364, 6.236000145419362, 6.457820872279776, 6.50070976857286, 6.961509673897431, 4.512213272296408, 4.0470174634855045, 2.432787270722598, 0.7713483310454722, 0.027332106079670626, 0.0]),
array([15.97064940217278, 16.86118082596039, 17.083453145161926, 14.117029932883897, 15.577215749027708, 14.10777642137027, 15.237678756707755, 15.801731508193013, 14.003039109833862, 15.843790536851358, 14.503616810005251, 14.192575051478379, 14.01318284967827, 14.764255645929136, 15.49950030291673, 14.2926495261976, 14.617632856708726]),
array([15.249604983379026, 15.893736251029765, 14.676663692876065, 14.089049634796357, 10.61240111889992]),
array([17.32904048145495]),
array([16.5827132679731, 17.220697874458228, 16.996277917183722, 16.012132631948685, 16.998142411698787, 16.997177606390242, 15.416489369148879, 14.782255407616315, 14.594549153089162, 15.923581283563276, 15.740219365147889, 14.450695101146435, 15.76018677531882, 14.795098252335162, 13.427600452972115, 12.931936732258334, 12.276725968039813, 11.440545262922553, 7.91172404311866, 10.515896295582515, 8.542060027574351, 8.850689872038618, 8.84442963832106, 10.328762865196062, 8.027673495088536, 11.314405302332496, 7.883601575292512, 9.933159163346158, 9.596375739011934, 9.67574212796524, 7.601243295726484, 8.122794708398104, 11.293994425068139, 8.925174288120013, 5.382547833630641, 7.105205138274203, 6.675600488638454, 7.173719668071294, 5.5403805855000785, 6.628676991497324, 5.4451108316969234, 5.568400235822478, 6.66556280272874, 6.5320074309030085, 5.957508711523804, 6.509925689023334, 5.977078313011514, 4.506675925974438, 3.946562338669114, 4.23449814731384, 2.6802122636304953, 2.898532391024133, 2.3653449854481927, 2.001371633287251, 0.7039589372892089, 0.0]),
array([16.36818126034146, 16.221724719739488, 16.087379366798007, 16.053554287894247, 16.69800943303938, 16.577779384696214, 16.556730233244775, 14.890890659654353, 14.678720403550253, 14.236100206297838, 15.129456363729672, 14.264073691289703, 15.46489134706612, 15.722639704820638, 14.829516665382721, 14.715239504815006, 14.259154415317084, 15.488554715585854, 15.269027267860803, 15.77056451177887, 13.859109326937796, 15.843683042082153, 14.879666402207153, 13.424460111958227, 13.771627459625671, 12.216897561161083, 10.11232411222147, 11.081449327033269, 9.480719971803362, 9.70478558242213, 10.118555981421899, 11.326016759536268, 9.023599316911277, 10.747229600451304, 10.602629881382885]),
array([16.618125252799597]),
array([16.10477880580164, 14.04444830563204, 13.903163226302054, 14.738339153320343, 13.89948946999278, 11.882020682680848, 13.475750665563616, 12.62059185536173, 7.83215724999498, 10.429278220107685, 10.378026086767486, 10.260282358690537, 7.490212134771643, 11.548847143048832]),
array([14.709783499646345]),
array([15.287466308351474, 15.228824716302196, 15.069437210506413, 14.332395755173494, 15.324323924122789, 14.457548695443835, 14.616466656820403, 15.309962037584995, 14.033934120255221, 14.770477308719618, 11.953600083609894, 12.332662098004105, 12.987936889355268]),
array([12.270188540593367]),
array([13.282199098309837, 7.280927918669421, 7.143907732061533, 5.5587883077642095, 4.318408950459682, 0.0]),
array([14.442255280102065, 14.304692895946289, 5.824270570639697, 5.577599608958032, 6.742543584828952, 4.535375578775537, 1.7347749078532009, 0.6958721322795081, 0.0]),
array([14.530819100715172, 14.10003598812218, 10.261070759664252, 6.106892590436804, 6.262701080334902, 6.724344434926037, 5.431170526388552, 6.295603645198284, 3.8953137991245237, 4.129268768754266]),
array([14.276278459545875, 13.83720597715741, 14.619964035735897, 14.463096875862863, 14.223316178782248, 14.020081981508607, 14.341137861169406, 14.257661706274757]),
array([12.106029549475185, 11.404463697003143, 11.202157202247802, 10.15658119541936]),
array([14.10391938929699, 13.86340377431577, 12.869441877306304, 10.431889315646131, 6.916022633182216, 6.572409716798435, 6.2558562671099285, 5.882943345139268, 4.607263633473453, 3.985768693296516, 2.8137486843729516, 1.9427164585776808, 0.3204922498083971, 0.09151468665846618, 0.0]),
array([13.967993557713902, 13.944628084217658, 14.042995285361215, 14.059237599099221, 13.935581434397465, 13.49334404415482, 13.786727909406538, 11.826028641710995, 12.833975015318243, 7.755541613607244, 7.970250444420207, 9.45788763806964, 11.220234923873214, 8.771789444754994, 10.852173481432493, 8.9982355542633, 9.453268618674999, 11.276250627256065, 11.402602350827143, 8.30802329339216, 8.202630744501938, 9.379366745617103, 8.484405225479806, 10.041370952569636, 11.360894988671523, 9.966648734314303, 10.80889397874304, 10.71250673091972, 8.94945273156881, 7.267135299506309, 8.601943668443425, 9.573886808406348, 5.557165582539058, 7.115233288412766, 7.134650929565739, 7.128271957696043, 6.366724717382101, 5.670204393583966, 6.541971757547096, 6.395767085553373, 5.454040141279187, 6.4749265332896355, 6.843586907139405, 6.461105678566064, 5.931568233420557, 6.856814853775382, 6.921108019556498, 6.810770854084672, 6.169605630044458, 3.672139301590954, 4.703811061932303, 5.216610146237098, 3.8297338672018717, 3.618387193933704, 4.012078279539276, 2.75898476818994, 3.0077345199523746, 2.40143571362921, 1.8775989526918582, 2.579641643089433, 1.9579692128919222, 1.9634819600666196]),
array([13.05629109452954, 8.621281104869466, 8.925581274755574, 7.413104027618399, 6.364996741825687, 2.3489466924709066]),
array([12.834747944702752]),
array([13.399019910007745, 12.622841783839513, 11.310118401648404, 8.894730366138106, 11.019701239626249, 5.732191962661818, 5.679328342731884, 5.486602586324173, 4.699453630479579, 4.536762720477908, 2.603425592382016, 3.2143384768107426, 1.853569952781728, 2.309769872971316, 2.461081032845034, 2.0260043892979924, 0.03441213963328586, 0.0]),
array([11.998786807969438]),
array([13.638113291867077, 11.973149909050331, 10.587860772136336, 11.188850827642069]),
array([13.320116683610442, 13.500728594424476, 13.50651007908477, 12.908702660535738, 10.486508103848518, 8.388738718563973, 8.216868161102882, 10.236147012512385, 7.7871718340041305, 5.763594313741626, 6.391129866836481, 4.163371413604397]),
array([12.261235873480617, 8.58931432234954, 10.593231458305155, 8.995077081732557, 5.777925287990457]),
array([12.681519272092949, 9.387946630103828, 9.984455675108917, 10.975701690262271]),
array([12.346072585906018, 12.92618502256251, 11.649242128388753, 11.234293361125127, 10.556648303596132, 10.09906767391389, 9.189134314455492, 8.37796540243246, 11.136484278309586, 11.576480070335464, 7.312815292655099, 6.639291244548282, 6.018289105073987, 7.1790923921254475, 7.039808186301356, 5.773497683692595, 5.4902359662898, 5.5475728404166125, 6.517105320531049, 5.574867091862346, 6.791693334612483, 4.417638873490139, 1.966398919239885, 2.358091130103327, 1.5438432756935874, 0.5234254799825109, 0.0]),
array([9.112053149394338, 10.512174904213591, 10.72364951993867, 10.681407281618691, 9.360061877347245, 11.225810921098022]),
array([12.783487656337128, 11.753679300482222, 13.029037023459606, 10.70351178524052, 10.914544492690354, 10.957067543057116]),
array([7.919725412951075, 10.08942869893573, 10.940503695163493, 9.328689401438984, 7.501538897845016, 8.361642438144864, 10.98695070286138, 8.297929868942678, 8.12067870707586, 10.962928309262136, 8.946011494858977, 10.535977199250956, 5.358014532422768, 6.540058929403791, 5.653702779874612, 6.669595599739085, 5.3649075506599555, 5.866026234725268, 5.606035478522283, 6.143023942816842, 3.984323109146571, 4.978154164158624, 4.09983336857097, 4.257409393364733, 4.576471711252576, 4.1830316979923445, 2.650437267277179, 2.976463258881579, 1.9540306071851496, 1.9557741153144916, 2.5093683948735905, 2.4064499084727125, 2.12847770567739, 1.8969979286716034, 0.2833718165512639, 0.5648182044311754, 0.417891555433155, 0.006538303003824311, 0.0]),
array([12.715534312968018, 12.517492510474598, 12.113303346377267, 11.234842928090728]),
array([12.35924121328646, 12.232656640458142]),
array([10.765226510899907, 10.38430106478448, 6.036923018135847, 6.040640635967703, 6.242523438666849, 5.866795461113823, 5.640120883555283, 5.634914347837455, 5.787265694727277, 4.917406096308203, 5.031120265327127, 3.2279241834851593, 3.4966084447617773, 3.445984763294497]),
array([9.507627865536413, 9.235035535861858, 8.42550380304379, 8.143943689642143, 6.337584957858715, 6.012492338172102, 6.021533520416464, 4.698462622906393, 3.352422535232964, 3.5812878961668737, 2.477601040267288, 2.1472951216379363, 1.48963616091755, 0.12290745290104378, 0.0]),
array([12.653944835796267, 12.489166226278494]),
array([8.209980948081641, 5.896392390373616, 5.44499240532379, 6.150725627370319, 5.124672469955403]),
array([12.049142206730089]),
array([11.809431949601874, 10.751609807558083, 7.870755056676883, 9.058442680483061, 9.569543921962422, 10.809884364897473, 11.005979697263665, 7.023642162513492, 7.022704522927887, 6.8636573374159315, 5.597279004702277, 5.603974441304587, 4.062004645509715]),
array([8.64639078448348, 7.880421553938826, 7.852056878232612, 6.370104806839846, 6.395182948516461, 6.303025601705379, 4.279419121900437, 3.0133176184734225, 3.4539525715027986, 2.326049775494784, 0.8788315932153459, 1.6103977205982218, 0.7428414005857494, 0.26447195526560585, 0.5276359790953868, 0.03532634406693788, 0.09823675253303812, 0.05871778075737914, 0.0]),
array([11.678658575770761, 10.32728598173543, 8.255572497467346, 8.597658845966329, 9.551404124819552, 7.740178694974528, 10.184058061905631, 9.899019829585638, 7.804915716201926, 9.965748023025853, 10.093989844501653, 8.555143999549927, 8.262057818353068, 9.76304742314043, 7.020223475709103, 7.115150587991732, 7.106057110623179, 7.14833543735888, 7.230340956723887]),
array([11.68897050190799, 10.727390253222369, 9.555668371976537, 9.024824798165652, 10.471710560812873, 10.556591103353144, 8.181730331229286, 7.248085968815677, 8.487865473938374, 11.241144136103536, 9.69454551074195, 10.566946175255616, 10.336466118660908, 8.255426043491433, 9.77075965141975, 8.705502336567472, 8.857742150444864, 11.052231756287345, 5.95901127455007, 7.0431127968552625, 6.326553892561134, 5.4469421099848585, 5.658798690537029, 6.332515730720683, 7.137145055785613, 5.728302284350516, 5.9055477238598115, 6.549617275879162, 6.397035596393443, 5.522923423811292, 6.922622803860578, 5.988811570623306, 5.85852266236086, 4.734385445690322, 4.749896232390599, 4.4267673652112975, 2.5926414296379683, 3.4389456446019486, 1.9936910036619362, 2.2729834286874873, 2.3310598314457605, 2.0772929873848103, 1.6937917191893168, 1.785367023203302, 0.12699375270036728, 0.09214913558924342, 0.0]),
array([5.627629627762782, 5.456365494732395, 6.930034424894617]),
array([8.28592044523757, 10.998930302299012, 8.846504052902786, 8.243224062489814, 9.335097721442642, 5.539381587789853, 5.493432359842545]),
array([9.780036090908167, 8.82934742559479, 11.009162215912792, 7.499206263772676, 7.320521423953777]),
array([10.729376199255457, 10.759826799082072, 10.768272566631078, 11.15996876575105]),
array([8.306714332643963, 9.52430468733147, 10.722870601450927, 9.253880178613501, 8.407363789686165, 9.810730825422326, 10.038952274858751, 7.997711831460183, 8.362772533695404, 8.33826006328519, 10.790277086834944, 9.948980450521573, 9.300334977508372, 8.31652862412934, 10.175834552029174, 8.801077492423195, 6.312139985186658, 5.5856807553378545, 6.452100546211021, 6.197843160488956, 6.379424736938726, 5.761458503430114, 6.598222929254773, 6.5500428230877725, 5.588555595891556, 5.490801178196633, 6.972336016318165, 5.37880572540294, 6.216414117331583, 6.757004949531645, 6.948547997411688, 5.845660003098827, 7.077282843744086, 5.9120583810125, 6.556041326726499, 6.852691395903936, 6.9188719795192934, 7.036372385648335, 6.277808805974882, 6.7225942529250835, 3.6274578162477358, 4.865279139620073, 5.168541536712025, 4.39195814800307, 4.619910515316686, 3.1911340544677995, 3.2152424857445543, 2.6683357370749476, 3.1968051640648425, 3.0632915591024124, 2.4410630338599026, 2.1804601138617974, 1.9464861955314547, 1.8748533502886553, 2.072079077836097, 2.060673585659349, 2.467904133431165, 2.1958456530800037, 1.582526598178212, 1.1135880779640868, 1.529676371542019, 1.711216696347992, 1.3823253113594012, 0.35942884885139553, 0.7499671533228757, 0.1558569244384188, 0.0]),
array([10.030770183731487]),
array([10.409652238252445, 9.887466632990108]),
array([8.957825476271244, 10.39611815157449, 7.428318243153607, 6.283523360777014, 5.7655401150677665, 6.607578699832991, 2.0383719312026254, 0.0]),
array([8.015795086397802, 7.523541930275658, 10.367152674791031, 8.944589396408627, 9.939135748821668, 6.85871259047355, 6.382089710989485]),
array([9.805636785825822, 7.76674249345754, 7.765800096274912, 9.194886173344996, 7.098666429328775, 7.018144354027143, 7.04499709017107, 6.74029040377751, 5.6877746156155995, 5.7848544070704575, 7.06404353521159, 4.74235926378196, 4.315349910687047, 4.0361352354828295, 3.945260012216231, 3.4006855405618444, 1.7530233297945046, 0.37349347526143767, 0.4746930510968856, 0.1027938860122425, 0.0]),
array([9.327923428014136]),
array([7.283515601897458, 9.250854035433433, 9.86092896275029, 9.453182779357858, 9.9969621637526, 8.248279691613044, 7.044039114474669, 6.058023124828932, 6.81653814360647, 7.235723289149516, 6.564241478796156, 5.799055026478649, 6.340547997834817, 6.170594152365258, 6.543409050882238, 6.289361795133132, 5.848505221032623, 6.749293737721213, 5.435515758149148, 6.205316530933343, 7.150609027645283, 5.739779938190032, 5.74536587198361, 6.477613967299157, 4.170031600084716, 5.230834363914172, 5.022876850626277, 5.257759326343481, 4.567791445651707, 4.544490313707034]),
array([8.759505335981189, 7.36461901871065, 7.1523386705475875]),
array([10.067350350186613]),
array([8.459833945008812, 7.5520456562875085, 6.5872327742933106, 6.742109884106176, 5.732921631568707, 6.418209351905604, 6.278884497081112, 6.646985542853872, 3.967411535355846, 4.872505529535077, 3.66074265838133, 2.9616253480898473, 2.701608586562744, 3.2612433039375732, 2.15599177814699, 1.8807817861831095, 2.3887411905522553, 2.0864819037586058, 0.5107087086497664, 0.0]),
array([9.53020863506264, 5.607504169898743, 2.375793411428291, 0.0]),
array([7.3725493621306235, 9.715594951408475, 7.594842107311393, 6.183932597262315, 6.324799815201966, 5.562523728212354, 7.038361407296382, 6.402690659899957, 5.805154858081664, 6.189813646231425, 6.8912392789987855, 5.941103635519457, 3.6279702078741636, 4.895977589176309, 4.282888949813445, 3.7480655981356725, 4.94592260701284, 3.958459440268208, 2.091809665314819, 2.3563941635480057, 2.3474309437641354, 2.347683882326324]),
array([9.15592407133658, 9.155973786917253, 9.569563805958195, 9.271821838434285, 5.7134675411811315, 5.539533911587886, 5.379832221922824, 6.439188503692528, 5.480231660813274, 6.833408354053946, 6.91368368785546, 6.64433137179064, 6.008731097600164, 5.879250132591078, 7.051306749892903, 5.356812570864752, 6.617458109550256, 6.6317265611288825, 6.636693989114947, 5.570633553823027, 5.771981216433253, 6.679235200035253, 3.637257169682635, 4.023386489498084, 5.2175005891531265, 4.086328709270041, 2.998320255904165, 3.4099424908497804, 2.846226589450856, 3.3274628307039773, 2.1515526269828116, 2.344077011216366, 2.1162322574166788, 1.8516336625726195, 2.4322392373416535, 2.0570167468590084, 2.5416820545129846, 2.4710279621656026, 2.3079277756808603, 1.3207089188092715, 0.9912065807076945, 1.2885129342818487, 0.29478260267210316, 0.015085848562996648, 0.0]),
array([8.498517961612396, 9.297319450819032, 9.39745829470992, 8.071918386264233, 5.985607768319149, 6.013357400122169, 6.642868025733913, 6.003026811963432, 5.7725270149711285, 5.810687837970053, 6.57399046996926, 6.038703627815946, 5.349068133397851, 6.108790550551253, 6.392229001206535, 6.493135956126645, 6.172279431291166, 5.420520765800315, 6.9027341839600584, 6.489475493542482, 3.846936432347789, 5.17102703810867, 4.759065081721364, 4.70055712102128, 4.964497065857015, 2.6669237264647743, 2.892775176820347, 3.3515910086058405, 2.027634003470333, 2.1528611792445846, 2.5787468425533713, 1.3396372133235976, 0.3528668435728274, 0.33069167377132413, 0.020466039860074975, 0.0]),
array([8.110682021029453, 3.3435981835196937, 2.353415543398962, 0.032643517831898936, 0.0]),
array([7.450318669951416, 8.91652948536748, 7.738885286840823, 8.238320864824013, 8.570425969822503, 8.649996214931432, 7.534465438842595, 8.107441239475172, 8.793583416663962, 7.750935232357897, 6.862744507988079, 5.621438840481159, 6.5450808727776355, 6.108015062460188, 5.5346981663786305, 5.65052831485084, 5.600892672705103, 6.86252556061042, 6.976593651193472, 6.960541380674129, 6.30904445172911, 5.333784065592536, 6.344882529500542, 5.974035050244735, 6.596146483812524, 5.825707716712455, 5.5235315492882995, 5.278644421947487, 4.255018592898834, 4.237159348344042, 4.665822633729921, 3.610040907837476, 4.9760819231174365, 3.1900341023363676, 2.6875948942895684, 1.8996054873155406, 2.300898838254461, 2.0856654173704725, 2.408984058897987]),
array([8.337605280912337, 7.86087347457081, 7.8356519860461695, 8.215394449233255, 7.750124903102366, 7.4515365362311226, 7.1178698549839]),
array([7.759675128442161, 8.135169096411992, 8.171546358697276, 8.340483250356169]),
array([5.7299016345150315, 6.739676353806444, 2.8998877522991755, 1.779228498605644, 0.3620701537707504, 0.057811281028043116, 0.0]),
array([7.716898262218723, 6.2148546949053936, 6.296326768140162, 6.630505259360083, 5.968064672811503, 7.195809756417701, 6.055324580593363, 5.550346438241898, 5.963147665555377, 3.392609453578069, 2.8773256560162426, 2.364094038569789, 2.057061362290328, 1.109293764960585, 0.06087729139549963, 0.1080669959647712, 0.1021810951437869, 0.01567439372914839, 0.0]),
array([6.558535421034792, 5.563117627060786, 6.626535739826762, 6.576688254945702, 6.202922304155184, 4.000293240112519, 2.5002378120253668, 1.916787241135732, 0.09287310511955647, 0.0]),
array([6.381225869719879, 5.374917980100767, 5.452742976840135, 6.025891870027235, 6.257372640858534, 3.3383087873555435, 1.0098023911284915, 0.8452139522385644, 0.4358944514547161, 0.0]),
array([7.781659956734707, 6.646658202695673, 6.437790335050774]),
array([6.725390700417272]),
array([7.869930953440461, 6.050408917584934, 5.939411169105757, 5.576386290995152, 7.073500493883538, 6.179904044133581, 6.967838629145659, 6.5672039738367936, 5.758973227834149, 5.192035527449128, 3.032932171005669]),
array([6.7620149423494285]),
array([7.022689107168308, 7.08326640224263, 6.534813179356269, 6.9589637323899485, 6.38102399214963]),
array([6.406219611554393, 7.0000997393935736, 5.81594559765318, 5.516236886498352, 6.7700517502758375, 4.697943695801359, 4.7003371120124475, 4.309221179644577, 2.7066990870859913, 3.036522521607511, 2.0621821865894905]),
array([7.307697547394198, 5.636901267090321, 5.753762044763293, 1.1733092856276823, 0.0]),
array([5.751996638071838, 4.090598447919315]),
array([6.588996490323819, 6.986672792796947, 6.782331198858488, 5.460483056269685, 6.471044685473719, 6.666157059181895, 5.845065950179947]),
array([6.608636880214746, 5.436605462247534, 6.926140729285028, 6.5791518044379975, 6.020235207410109, 3.7792468820504377, 4.928115611193421, 5.303836266342702, 2.213807426185794, 2.345247334491411, 2.302986268540491, 1.4920211593221506, 0.4983592674407002, 0.6172186819195468, 0.0]),
array([5.951731991032638, 0.11395148584514009, 0.0]),
array([6.740405655033249, 0.0]),
array([5.6210130616082745, 3.954751950220655, 4.201558433018315]),
array([6.574845783729701, 6.637855517414757, 4.407562185372881, 4.097603783490713, 2.8175494116653095, 2.658334599986281, 2.2734816739117516, 0.25694638309271944, 0.0]),
array([5.5187733326630655, 3.1679561018923557, 0.6142111875738092, 0.0]),
array([5.82502906135042, 5.657328880244274, 5.824724628906794, 6.549412180079368, 6.286281791271158]),
array([5.963260222900441]),
array([6.071022147383441, 5.695916703438825, 2.2504315816360316]),
array([5.654855192018655, 5.696317117519024, 5.763453841353742, 5.472800999251581, 5.673661140573408, 4.114280925515776, 5.093267126935917, 4.452567775881951]),
array([5.606535860922437]),
array([4.62886190815635, 2.2207488481549413, 1.136010795388155, 0.0]),
array([5.480549134045719, 5.470644487737925, 4.869439493763893, 4.82404675763039, 4.830409731634978, 3.089506179853174, 3.386111343143156, 0.37957438725275694, 0.0]),
array([5.4011039191689045, 5.6419765405689635, 5.619952083043824, 5.501945674971706, 4.068549471393507, 3.6423616794659495, 2.6958022400600363, 3.0748307484543336, 1.9431971436182267, 2.2570654504015786, 0.359673117322403, 0.1354485851166498, 0.17231780049789647, 0.04091298124645787, 0.03953904345143895, 0.04692205000918884, 0.01309289067000613, 0.0]),
array([5.3350385284416335, 5.686859447500556, 5.628557036116136, 2.707163728052573, 3.0715526726556437, 2.6330625387771827, 2.3065953026762602, 1.9906853219183873, 1.2463475703064024, 1.6951211296854423, 1.7080155702536057, 0.3270420560548004, 0.09603074152900258, 0.0]),
array([5.539329934061675, 4.69284274959425]),
array([5.586063412540089, 5.200154018869868, 5.13449943163119, 3.645305720033886, 2.6088613580783457, 1.9442209956561243, 2.0186297871290027, 1.9005695133577354, 2.0251866142616, 0.6312517197623393, 0.4082156218222346, 0.6755373204975245, 0.0]),
array([3.019821483191276, 2.382720979487467, 2.0529480042150183, 1.3582672876210697, 0.0]),
array([5.1361018476733955]),
array([5.016431048519445, 3.6813635529789286, 5.019186850111188, 4.167386915596892, 4.793304778131502, 2.8454414205486147, 1.8538115758196856, 1.841346650101995, 2.269016549458293, 2.0226587278429067, 1.8870901160278213, 2.5313552536160695, 0.0626372676536895, 0.0]),
array([0.2636461598642573, 0.0]),
array([4.206089228834379, 4.454528870999957, 4.326739481703524, 2.5169386221519026, 2.4478722758223084, 0.3514065760581768, 0.0]),
array([3.898150695813649, 2.050143903918714, 0.05627913180306561, 0.0]),
array([4.0701235470487305, 4.364321043666119, 1.8440201186543839, 0.09229200247516173, 0.013631216525692713, 0.0]),
array([2.9602538451962026, 0.6654634829728329, 0.1654868942421811, 0.0]),
array([2.3860381293566593, 2.166561172916067, 2.0219546837273796, 0.7981258562463622, 0.5845936640589982, 0.0]),
array([0.11479883101910175, 0.0]),
array([3.979938947669553, 3.4635667765146945, 2.1160816116945527, 2.205904334995202, 0.06337125013866161, 0.0]),
array([3.3923698713972494, 1.826321730122811, 2.24276538297535, 2.0307232316517423, 2.1357349579084546, 1.9380668812372952, 1.6971722309426578, 0.5506380783868909, 0.39044470515003943, 0.0]),
array([2.3949971114676236]),
array([2.7762363855715244, 1.1397819026382745, 0.0]),
array([3.7341484429751266, 3.367596201582099, 3.413216038168063, 3.4658329922220275, 2.032558489005292, 2.2016870972273574, 1.9647443579141748]),
array([1.1670913858695584, 0.6861294356075638, 0.0]),
array([3.3514247478362202, 2.7949294578936534, 2.974014317791399, 3.4805617652990857, 0.38021774003297604, 0.13370477129007952, 0.0]),
array([3.7308950207539193]),
array([3.497420998936531, 3.018582804535545, 1.837375320368736, 2.5759149610772205, 0.41750599471112315, 0.06989833171084917, 0.0]),
array([3.206932928183045, 3.483616435763804, 1.9614209962233313, 1.8557389679367695, 1.245644221791263, 1.6778250751542816, 0.7405757598544719, 0.0]),
array([3.240773840138124, 1.9824927388957607, 0.2603100462795557, 0.0]),
array([2.256087550292926, 2.3030974975769127, 0.9635020682792323, 0.4213968412714198, 0.0830977653691108, 0.0]),
array([2.883047892518001, 2.9856424741778427, 2.8250890461249174, 2.5806302462994406, 2.319367149515724, 2.0128210399349302, 2.3275098094424607, 2.1040930761500327, 1.3170944551056079, 1.5533397114497989, 0.8874522895799165, 0.1452447018671721, 0.3115484768251454, 0.0]),
array([2.092897137826998, 0.0]),
array([2.587595201623369]),
array([2.550678354401193, 0.9512570739252721, 0.17803748393252772, 0.0]),
array([2.3128090556492005]),
array([0.7494158429165156]),
array([2.4058322863998898, 2.201154250199644, 1.9910894093133198, 2.4113981574649506, 0.060090825714591664, 0.0]),
array([2.534601211812673, 0.1525727269929239, 0.6418309256491892, 0.062042856958573006, 0.0]),
array([2.378856846602018, 0.37363627060663523, 0.3266226433824036, 0.6362062987538184, 0.07590297871742135, 0.0]),
array([2.2849565338288502, 2.1447040012497487, 0.9418086310011653, 0.0]),
array([0.31606038680629184, 0.0]),
array([2.0388895818458193, 0.0]),
array([2.1871644879443592, 0.8284693733118803, 0.6930008783299166, 0.0]),
array([1.9088809388229337, 0.6749787400500391]),
array([2.2004978826465664, 2.1853739246741886, 0.994176887272215, 1.2335564528054521, 0.1671432784370075, 0.5640048671819172, 0.4045744925242921, 0.12418610488636388, 0.0]),
array([1.411303532739621, 0.10643958983644514, 0.06364237582482116, 0.0]),
array([0.48866518790839164, 0.7546978885782382, 0.7732558212774939, 0.2511773122771406, 0.06989807425357161, 0.09481027252748798, 0.0]),
array([1.3700484239024084, 0.0]),
array([1.5732377513875895, 0.0]),
array([1.1759714984692953, 0.12153492072229367, 0.0]),
array([1.5517517192491932, 0.6222954733568804, 0.0]),
array([1.401556277211617, 1.3829356699961146, 0.0]),
array([0.1859690808311908, 0.0]),
array([0.03734453967591343, 0.0]),
array([0.6038316461118178, 0.0]),
array([0.21832355025148165, 0.1835306133283091, 0.0]),
array([0.8043544212358469, 0.5324226547175702, 0.0]),
array([0.4823929014475356, 0.0]),
array([0.5557414287367697, 0.0]),
array([0.34612138707478896, 0.0]),
array([0.10042684549180035, 0.0]),
array([0.5443876446172776, 0.0]),
array([0.6007550525201222, 0.0]),
array([0.3152332121930896, 0.6995767834903662, 0.23436563906019947, 0.5333952523317389, 0.0]),
array([0.42276775899920577, 0.00330078197227536, 0.0]),
array([0.8047670337369227, 0.06395857937651611, 0.09782292100852893, 0.11407126691366809, 0.0]),
array([0.43726759310770486, 0.0]),
array([0.6786099415421072, 0.08009543029681783, 0.0]),
array([0.5807687106527918, 0.0]),
array([0.1871438063500392, 0.0]),
array([0.07114587335841972, 0.08446935458516729, 0.0]),
array([0.1393817331717854, 0.0]),
array([0.29422615500787513, 0.07032507132408633, 0.0]),
array([0.007194585116282404, 0.09828323324196382, 0.0]),
array([0.14654499853382708, 0.0]),
array([0.151698156845015, 0.0]),
array([0.07648546513053375, 0.0])
]
d = [data_1]
names = ["83"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T48', 'T49', 'T50', 'T53', 'T54', 'T55', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T68', 'T69', 'T70', 'T71', 'T72', 'T74', 'T75', 'T76', 'T78', 'T79', 'T81', 'T83', 'T84', 'T85', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T93', 'T95', 'T96', 'T97', 'T98', 'T100', 'T101', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T131', 'T132', 'T134', 'T135', 'T136', 'T137', 'T138', 'T139', 'T140', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T149', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T161', 'T162', 'T163', 'T164', 'T165', 'T167', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T180', 'T182', 'T183', 'T185', 'T186', 'T187', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T198', 'T200', 'T201', 'T204', 'T205', 'T206', 'T208', 'T209', 'T210', 'T211', 'T213', 'T215', 'T218', 'T220', 'T222', 'T223', 'T224', 'T226']
def get_taxa_names(): return taxa_names