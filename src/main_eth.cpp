
#include <TP.hpp>
#include <Util.hpp>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/matrix_decompose.hpp>

#include <format>
#include <charconv>
#include <future>

#include <raylib.h>
#include <raymath.h>

#include <imgui.h>
#include <imgui_impl_raylib.h>

#include <font_regular.h>
#include <font_bold.h>
#include <logo.h>

int main() {

    using namespace TP;

    int32_t wind_w = 800;
    int32_t wind_h = 600;
    if(std::filesystem::exists("config.json")){
        std::ifstream i("config.json");
        json o;
        i >> o;
        wind_w = o["width"].template get<int32_t>();
        wind_h = o["height"].template get<int32_t>();
    }

    SetConfigFlags(FLAG_WINDOW_RESIZABLE); 
    SetConfigFlags(FLAG_MSAA_4X_HINT); 
    InitWindow(wind_w, wind_h, "TurboPacker");

    Image logo;
    logo.data = (void*)logo_data;
    logo.width = 301;
    logo.height = 335;
    logo.mipmaps = 1;
    logo.format = PixelFormat::PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;

    SetWindowIcon(logo);  

    const auto scale = GetWindowScaleDPI();

    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; 
    ImFont* regular = io.Fonts->AddFontFromMemoryTTF(Montserrat_Regular_ttf, sizeof(Montserrat_Regular_ttf), 14.f * std::max(scale.x, scale.y));
    ImFont* regular_l = io.Fonts->AddFontFromMemoryTTF(Montserrat_Regular_ttf, sizeof(Montserrat_Regular_ttf), 16.f * std::max(scale.x, scale.y));
    ImFont* bold = io.Fonts->AddFontFromMemoryTTF(Montserrat_Bold_ttf, sizeof(Montserrat_Bold_ttf), 16.f * std::max(scale.x, scale.y));
    ImFont* massive = io.Fonts->AddFontFromMemoryTTF(Montserrat_Bold_ttf, sizeof(Montserrat_Bold_ttf), 20.f * std::max(scale.x, scale.y));
    ImGui::StyleColorsDark();
    
    ImGui_ImplRaylib_Init();
    Imgui_ImplRaylib_BuildFontAtlas();

    auto conf = Disk::load<float, CostFunction::CF_Krass, uint16_t, uint32_t, 15, std::allocator<uint16_t>>(std::filesystem::path("config.json"));
    conf.BoxType = TP::Detail::BoxGenerationType::LIST;
    conf.LookAheadSize = 50; 
    conf.EmptryTries = 500;
    conf.MaxEmptryTries = 1000;

    conf.Bins.clear();
    //E
    conf.Bins.push_back({ {75.f, 50.f}, 200 }); //0
    conf.Bins.push_back({ {75.f, 150.f}, 50 }); //1
    conf.Bins.push_back({ {75.f, 150.f}, 50 }); //2
    conf.Bins.push_back({ {75.f, 200.f}, 50 }); //3

    //T
    conf.Bins.push_back({ {75.f, 50.f}, 200 }); //4
    conf.Bins.push_back({ {75.f, 200.f}, 50 }); //5

    //H
    conf.Bins.push_back({ {75.f, 50.f}, 250 }); //6
    conf.Bins.push_back({ {75.f, 125.f}, 50 }); //7
    conf.Bins.push_back({ {75.f, 50.f}, 250 }); //8

    using Prmise = decltype(solve(conf));
    std::unique_ptr<Prmise> pr = nullptr;

    //-------------------------------------
    Util::Camera3d camera;
    camera.up = { 0., 1., 0. };
    camera.target = { 0.f, 0.f, 0.f };
    camera.camDist = 800.f;
    camera.tiltAngle = -65.f;
    //-------------------------------------

    std::vector<TP::Detail::Entry<float>> b;
    int32_t cc = 0;

    const auto col = [&](float _scalar) -> Color {
        const float RedSclr = std::clamp<float>((1.0f - _scalar)/0.5f,0.f,1.f);
        const float GreenSclr = std::clamp<float>((_scalar/0.5f),0.f,1.f);
        const uint8_t R = (uint8_t)(255 * RedSclr);
        const uint8_t G = (uint8_t)(255 * GreenSclr);
        const uint8_t B = 0;
        return Color{ R, G, B, 175 };
    };

    const auto col_img = [&](float _scalar) -> ImVec4 {
        const auto cl = col(_scalar);
        constexpr float scl = 1.f / 255.f;
        return ImVec4{ float(cl.r) * scl, float(cl.g) * scl, float(cl.b) * scl, 1.f };
    };

    const char* str_modes[] = { "Random", "List" }; //, "Validate" };
    int cur_mode = (int)conf.BoxType;
    uint64_t seedval = conf.Seed;
    std::unique_ptr<Util::Detail::EvalRes> eval_res = nullptr;

    //-------------------------------------

    while (!WindowShouldClose()) {

        if(IsWindowResized()){
            wind_w = GetScreenWidth(), 
            wind_h = GetScreenHeight();
        }

        ImGui_ImplRaylib_ProcessEvents();
       
        //------------------------

        ImGui_ImplRaylib_NewFrame();
        ImGui::NewFrame();
        ImGui::Begin("PackerWidget", (bool*)nullptr, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize);
        //-------------------------------------

        {//Time label
            const auto txt = std::format("{}s", pr ? pr->getTime() : 0.f); 
            const ImVec2 textSize = ImGui::CalcTextSize(txt.data());
            const float centeredX = (ImGui::GetContentRegionAvail().x - textSize.x) * 0.5f;
            ImGui::PushFont(massive);
            if(!pr || pr->isDone()) ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 1.0f, 0.0f, 1.0f));
            else ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
            ImGui::SetCursorPosX(centeredX);
            ImGui::Text(txt.data());
            ImGui::PopStyleColor(1);
            ImGui::PopFont();
            ImGui::Dummy({0, 10});
            ImGui::Separator();
            ImGui::Dummy({0, 10});
        }

        {//progress
            ImGui::PushFont(regular_l);
            const ImVec2 ts1 = ImGui::CalcTextSize("Total1000");
            ImGui::BeginChild("TotalMain", ts1, false);
            ImGui::Text("Total");
            ImGui::EndChild();
            ImGui::SameLine();   
            const ImVec2 ts2 = ImGui::CalcTextSize("1000");
            ImGui::BeginChild("ProgressMain", ts2, false);
            ImGui::Text(std::format("{}", pr ? pr->getTotalBoxCount() : 0).data());
            ImGui::EndChild();
            ImGui::SameLine();
            const float prg = pr ? pr->getTotalPackDensity() : 0.f;
            ImGui::PushStyleColor(ImGuiCol_PlotHistogram, col_img(prg));
            ImGui::ProgressBar(prg);
            ImGui::PopStyleColor(1);
            ImGui::PopFont();
            ImGui::Dummy({0, 6}); 

            //bin progress
            const ImVec2 ts3 = ImGui::CalcTextSize("1000");
            for(size_t i = 0; i < conf.Bins.size(); ++i) {
                ImGui::BeginChild(std::format("Total##{}", i).data(), { ts1.x, ts3.y }, false);
                ImGui::Text(std::format("Bin {}", i).data());
                ImGui::EndChild();
                ImGui::SameLine();
                const ImVec2 ts4 = ImGui::CalcTextSize("1000");
                ImGui::BeginChild(std::format("Progress##{}", i).data(), { ts2.x, ts4.y }, false);
                ImGui::Text(std::format("{}", pr ? pr->getBoxCount(i) : 0).data());
                ImGui::EndChild();
                ImGui::SameLine();
                const float prg = pr ? pr->getPackDensity(i) : 0.f;
                ImGui::PushStyleColor(ImGuiCol_PlotHistogram, col_img(prg));
                ImGui::ProgressBar(prg);
                ImGui::PopStyleColor(1);
            }
        }

        ImGui::Dummy({0, 10});
        ImGui::Separator();
        ImGui::Dummy({0, 10});

        //-------------------------------------
        if(!pr || pr->isDone()){
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.f, 1.f, 0.f, 1.f));
            if(ImGui::Button("Pack", ImVec2(100, 30))) {

                if(conf.BoxType == Detail::BoxGenerationType::VALIDATE && eval_res == nullptr) {
                    std::cerr << std::format("Generate a solution first!") << std::endl;
                } else if(conf.BoxType == Detail::BoxGenerationType::VALIDATE && eval_res != nullptr) {
                    Config<float, CostFunction::CF_Krass> cc = conf;
                    cc.BoxType = Detail::BoxGenerationType::LIST;

                    std::sort(eval_res->decomp_.begin(), eval_res->decomp_.end(), [](const auto& _v1, const auto& _v2){
                        return _v1.first < _v2.first;
                    });

                    Detail::BoxList<float> list;
                    list.shuffle_boxes_ = true;
                    //for(const auto& [id, bx] : eval_res->decomp_){
                    //    list.list_.push_back( { {bx.getSize().x, conf.Height - 2.f, bx.getSize().y }, 1, TP::PF_Y_XZ | TP::PF_Y_ZX });
                    //}
                    cc.BoxList = list;

                    pr = std::make_unique<Prmise>(solve(cc));
                } else pr = std::make_unique<Prmise>(solve(conf));
                
            }
            ImGui::PopStyleColor(2);

            if(!b.empty()){
                ImGui::SameLine();
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1.f, 0.f, 0.0f, 1.f));
                if(ImGui::Button("Reset", ImVec2(100, 30))) {
                    b.clear();
                    pr = nullptr;
                }
                ImGui::PopStyleColor(1); 
            } 

        } else {
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1.f, 0.f, 0.0f, 1.f));
            if(ImGui::Button("Stop", ImVec2(100, 30))) {
                pr->stop(); 
            }
            ImGui::PopStyleColor(1);
        }
        //-------------------------------------
        ImGui::Dummy({0, 10});
        ImGui::Separator();
        ImGui::Dummy({0, 4});
        ImGui::PushFont(bold);
        ImGui::Text("Config");
        ImGui::PopFont();
        ImGui::Dummy({0, 4});

        if(ImGui::BeginCombo("Mode", str_modes[cur_mode])) {
            for(int32_t i = 0; i < IM_ARRAYSIZE(str_modes); ++i){
                if(ImGui::Selectable(str_modes[i], cur_mode == i)){
                    cur_mode = i;
                    switch(cur_mode){
                        case 0: conf.BoxType = Detail::BoxGenerationType::RANDOM; break;
                        case 1: conf.BoxType = Detail::BoxGenerationType::LIST; break;
                        case 2: conf.BoxType = Detail::BoxGenerationType::VALIDATE; break;
                    }
                    b.clear();
                    pr = nullptr;
                } 
                if(i == cur_mode) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        ImGui::Checkbox("Multithreading", &conf.MultiThreading);
        if(conf.MultiThreading) ImGui::InputInt("Threads", (int32_t*)&conf.NumThreads);
        ImGui::Checkbox("Random Seed", &conf.UseRandomSeed);
        {
            seedval = conf.Seed;
            char buf[20];
            std::sprintf(buf, "%llu", seedval);
            if(ImGui::InputText("Seed", buf, IM_ARRAYSIZE(buf), ImGuiInputTextFlags_CallbackAlways, [](ImGuiInputTextCallbackData* data) {
                if (data->EventFlag == ImGuiInputTextFlags_CallbackAlways) {
                    uint64_t* p_value = (uint64_t*)data->UserData;
                    std::string input_text(data->Buf);
                    uint64_t value;
                    auto [ptr, ec] = std::from_chars(input_text.data(), input_text.data() + input_text.size(), value);
                    if (ec != std::errc::invalid_argument && ec != std::errc::result_out_of_range)
                        *p_value = value; 
                }
                return 0;
            }, &seedval)){
                if(!conf.UseRandomSeed) conf.Seed = seedval;
            }
        }
        ImGui::Checkbox("Allow Overlap", &conf.AllowOverlap);

        ImGui::Dummy({0, 4});
        ImGui::Separator();
        ImGui::Dummy({0, 4});

        for(size_t i = 0; i < conf.Bins.size(); ++i){
            ImGui::Text(std::format("Bin {}", i).data());
            ImGui::InputFloat2(std::format("Bounds##{}", i).data(), glm::value_ptr(conf.Bins[i].Bounds));
            ImGui::InputInt(std::format("Height##{}", i).data(), (int32_t*)&conf.Bins[i].Height);
            ImGui::Dummy({0, 3});
        }

        ImGui::Dummy({0, 4});
        ImGui::Separator();
        ImGui::Dummy({0, 4});

        ImGui::InputInt("Empty Tries", (int32_t*)&conf.EmptryTries);
        ImGui::InputInt("Max Empty Tries", (int32_t*)&conf.MaxEmptryTries);
        ImGui::InputInt("Look Ahead", (int32_t*)&conf.LookAheadSize);

        ImGui::Dummy({0, 4});
        ImGui::Separator();
        ImGui::Dummy({0, 4});

        if(conf.BoxType == Detail::BoxGenerationType::RANDOM){
            ImGui::Checkbox("Use Cubes", &conf.CubeRandomBoxes);
            ImGui::InputFloat("Min Box Volume", &conf.MinBoxVolume);
            ImGui::InputFloat("Max Box Volume", &conf.MaxBoxVolume);
        } else if(conf.BoxType == Detail::BoxGenerationType::LIST) {
            ImGui::Checkbox("Enforce Misses", &conf.EnforceMisses);
        } else if(conf.BoxType == Detail::BoxGenerationType::VALIDATE) {

            ImGui::InputInt("Box Count", (int32_t*)&conf.EvalBoxCount);
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.f, 1.f, 0.f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
            // if(ImGui::Button("Generate\nProblem", ImVec2(100, 60))) {
            //     if(conf.UseRandomSeed){
            //         conf.Seed = std::random_device()();
            //     }
            //     eval_res = Util::create_ground_truth(glm::vec<2, int32_t>{(int32_t)conf.Bounds.x, (int32_t)conf.Bounds.y}, conf.EvalBoxCount, conf.Seed);
            // }

            // ImGui::PopStyleColor(2);
            // ImGui::SameLine();
            // ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1.f, 0.f, 0.0f, 1.f));
            // if(ImGui::Button("Reset", ImVec2(100, 60))) {
            //     eval_res = nullptr;
            // }
            // ImGui::PopStyleColor();

            // if(eval_res){
            //     const float ww = (ImGui::GetContentRegionAvail().x - 200) * 0.5f;
            //     const float ar = conf.Bounds.y / conf.Bounds.x;
            //     ImGui::Dummy({0, 10});
            //     ImGui::SetCursorPosX(ww);
            //     ImGui::Image((void*)&eval_res->tex_, 
            //         ImVec2(200, int32_t(200.f * ar))
            //     );
            // }
            
        }

        if(!ImGui::IsWindowHovered()) camera.update(GetFrameTime());
       
        ImGui::End();
        ImGui::Render();
        
        //------------------------

        BeginDrawing();
        ClearBackground(WHITE);

        BeginMode3D(camera);
        //DrawGrid(100, 10.0f);

        for(size_t i = 0; i < 9; ++i){

            float dx = 0.f;
            float dy = 0.f;
            float dz = 0.f;

            switch(i) {
                case 0:
                    dx = 0.f;
                    dy = 0.f;
                    dz = 0.f;
                break;
                case 1:
                    dx = 0.f;
                    dy = 50.f;
                    dz = 0.f;
                break;
                case 2:
                    dx = 0.f;
                    dy = 50.f;
                    dz = 100.f;
                break;
                case 3:
                    dx = 0.f;
                    dy = 0.f;
                    dz = 200.f;
                break;
                case 4:
                    dx = 0.f;
                    dy = 325.f;
                    dz = 0.f;
                break;
                case 5:
                    dx = 0.f;
                    dy = 250.f;
                    dz = 200.f;
                break;
                case 6:
                    dx = 0.f;
                    dy = 525.f;
                    dz = 0.f;
                break;
                case 7:
                    dx = 0.f;
                    dy = 575.f;
                    dz = 100.f;
                break;
                case 8:
                    dx = 0.f;
                    dy = 700.f;
                    dz = 0.f;
                break;
            }

            const auto& bin = conf.Bins[i];

            const auto ext = bin.Bounds * 0.5f;
            DrawCubeWiresV(Vector3{ext.y + dy, bin.Height * 0.5f + dz, ext.x + dx}, Vector3{bin.Bounds.y, (float)bin.Height, bin.Bounds.x}, PURPLE);

            //if(i == 0){
            //    DrawLine3D(Vector3{0, 0, 0}, Vector3{(float)bin.Height, 0, 0}, RED);
             //   DrawLine3D(Vector3{0, 0, 0}, Vector3{0, (float)bin.Height, 0}, GREEN);
            //    DrawLine3D(Vector3{0, 0, 0}, Vector3{0, 0, (float)bin.Height}, BLUE);
            //}

        }

        if(pr && cc++%30 == 0) b = pr->data_cpy();  

        for (size_t i = 0; i < b.size(); ++i) {
            const auto& e = b[i];    

            float dx = 0.f;
            float dy = 0.f;
            float dz = 0.f;

            switch(e.bin_id_) {
                case 0:
                    dx = 0.f;
                    dy = 0.f;
                    dz = 0.f;
                break;
                case 1:
                    dx = 0.f;
                    dy = 50.f;
                    dz = 0.f;
                break;
                case 2:
                    dx = 0.f;
                    dy = 50.f;
                    dz = 100.f;
                break;
                case 3:
                    dx = 0.f;
                    dy = 0.f;
                    dz = 200.f;
                break;
                case 4:
                    dx = 0.f;
                    dy = 325.f;
                    dz = 0.f;
                break;
                case 5:
                    dx = 0.f;
                    dy = 250.f;
                    dz = 200.f;
                break;
                case 6:
                    dx = 0.f;
                    dy = 525.f;
                    dz = 0.f;
                break;
                case 7:
                    dx = 0.f;
                    dy = 575.f;
                    dz = 100.f;
                break;
                case 8:
                    dx = 0.f;
                    dy = 700.f;
                    dz = 0.f;
                break;
            }  

            const Vector3 ns = {glm::length(glm::vec3(e.tf_[0])), glm::length(glm::vec3(e.tf_[1])), glm::length(glm::vec3(e.tf_[2])) };
            const double temp = 1. - 1. / (conf.MaxBoxVolume - conf.MinBoxVolume) * ((ns.x * ns.y * ns.z) - conf.MinBoxVolume);
            const auto trans = glm::vec<3, float>(e.tf_[3]) + glm::vec<3, float>(dy, dx, dz);
            const Color cc = conf.BoxType == Detail::BoxGenerationType::VALIDATE && eval_res ? 
                Color{ eval_res->colmap_[e.id_].r, eval_res->colmap_[e.id_].g, eval_res->colmap_[e.id_].b, 175 } : col(temp);
            DrawCubeV(Vector3{ trans.x, trans.z, trans.y }, Vector3{ ns.x, ns.z, ns.y }, cc);
            DrawCubeWiresV(Vector3{ trans.x, trans.z, trans.y }, Vector3{ ns.x, ns.z, ns.y }, BLACK);
        }

        EndMode3D();
      
        ImGui_ImplRaylib_RenderDrawData(ImGui::GetDrawData());
        DrawFPS(GetScreenWidth() - 100, 5);
        EndDrawing();

    }
    ImGui_ImplRaylib_Shutdown();
    CloseWindow();

    json config = Disk::save(conf);
    config["width"] = wind_w;
    config["height"] = wind_h;
    std::ofstream s("config.json");
    s << std::setw(4) << config;
    return 0;
}
