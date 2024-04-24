#pragma once

#include <format>
#include <string>
#include <filesystem>
#include <fstream>

#include <raylib.h>
#include <raymath.h>

#include <nlohmann/json.hpp>

#include <TP.hpp>

namespace Util {

    namespace Detail {

        struct EvalRes {
            glm::vec<2, int32_t> bounds_;
            std::vector<Color> colmap_;
            std::vector<std::pair<int32_t, TP::Detail::FBox<2, int32_t>>> decomp_;
            Texture2D tex_;
            TP::Detail::BoxList<float> boxlist_;
        };//EvalRes

    }//Detail

    [[nodiscard]] inline std::unique_ptr<Detail::EvalRes> create_ground_truth(
        const glm::vec<2, int32_t> _bounds,
        const int32_t _num_boxes,
        const uint64_t _seed
    ) {

        std::mt19937_64 rng (_seed);

        std::unique_ptr<Detail::EvalRes> out = std::make_unique<Detail::EvalRes>();
        out->bounds_ = _bounds;

        out->colmap_.resize(_num_boxes);
        for(int32_t i = 0; i < _num_boxes; ++i)
            out->colmap_[i] = ColorFromHSV(i * (360.f / (_num_boxes + 1)), 1.f, 1.f);

        std::vector<std::pair<int32_t, TP::Detail::Interval<float>>> bxs;
        bxs.resize(_num_boxes);

        const float median = out->bounds_.x * out->bounds_.y / _num_boxes;
        for(int32_t i = 0; i < _num_boxes; ++i)
            bxs[i] = { i, {median * 0.5f, median * 1.5f} };

        out->decomp_ = TP::squarify(TP::Detail::FBox<2, float>{ {0.f, 0.f}, { float(_bounds.x), float(_bounds.y) } }, bxs, rng());

        std::vector<unsigned char> img_data;
        img_data.resize(4 * _bounds.x * _bounds.y, 0);

        for(const auto[idx, b] : out->decomp_){

            //fill
            for(size_t y = b.min_.y; y < b.max_.y; ++y){
                for(size_t x = b.min_.x; x < b.max_.x; ++x){
                    const size_t i = x + y * _bounds.x;
                    if(4*i >= img_data.size()) continue;
                    img_data[4*i] = out->colmap_[idx].r;
                    img_data[4*i + 1] = out->colmap_[idx].g;
                    img_data[4*i + 2] = out->colmap_[idx].b;
                    img_data[4*i + 3] = 255;
                }
            }

            //border y
            for(size_t y = b.min_.y; y < b.max_.y; ++y){
                const size_t i1 = b.min_.x + y * _bounds.x;
                const size_t i2 = (b.max_.x - 1) + y * _bounds.x;
                if(4*i1 >= img_data.size()) continue;
                if(4*i2 >= img_data.size()) continue;
                img_data[4*i1] = 0;
                img_data[4*i1 + 1] = 0;
                img_data[4*i1 + 2] = 0;
                img_data[4*i1 + 3] = 255;

                img_data[4*i2] = 0;
                img_data[4*i2 + 1] = 0;
                img_data[4*i2 + 2] = 0;
                img_data[4*i2 + 3] = 255;
            }

            //border x
            for(size_t x = b.min_.x; x < b.max_.x; ++x){
                const size_t i1 = x + b.min_.y * _bounds.x;
                const size_t i2 = x + (b.max_.y - 1) * _bounds.x;
                if(4*i1 >= img_data.size()) continue;
                if(4*i2 >= img_data.size()) continue;
                img_data[4*i1] = 0;
                img_data[4*i1 + 1] = 0;
                img_data[4*i1 + 2] = 0;
                img_data[4*i1 + 3] = 255;

                img_data[4*i2] = 0;
                img_data[4*i2 + 1] = 0;
                img_data[4*i2 + 2] = 0;
                img_data[4*i2 + 3] = 255;
            }

        }

        Image img;
        img.mipmaps = 0;
        img.width = _bounds.x;
        img.height = _bounds.y;
        img.format = PixelFormat::PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
        img.data = img_data.data();

        out->tex_ = LoadTextureFromImage(img);

        return out;
    }//create_groundtruth

    template<class T = float>
    class Camera3d : public Camera3D {

        Vector2 cursorPos = GetMousePosition();
        bool isDirty = true;

    public:

        T camDist = 250.; 
        T rotAngle = 45; 
        T tiltAngle = 4; 
        T rotSpeed = 0.5; 
        T moveSpeed = 50.; 
        T zoomSpeed = 15.;

        Camera3d() {
            projection = CAMERA_PERSPECTIVE;
            fovy = 65.f;
        }

        void update(const T _dt) {

            if (IsMouseButtonDown(1)) {
                Vector2 newPos = GetMousePosition();

                rotAngle -= (newPos.x - cursorPos.x) * rotSpeed;
                tiltAngle -= (newPos.y - cursorPos.y) * rotSpeed;

                if (tiltAngle > 89)
                    tiltAngle = 89;
                if (tiltAngle < -89)
                    tiltAngle = -89;

                isDirty = true;
            }

            cursorPos = GetMousePosition();

            Vector3 moveVec = { 0, 0, 0 };

            if (IsKeyDown(KEY_W))
                moveVec.z = -moveSpeed * _dt;
            if (IsKeyDown(KEY_S))
                moveVec.z = moveSpeed * _dt;

            if (IsKeyDown(KEY_A))
                moveVec.x = -moveSpeed * _dt;
            if (IsKeyDown(KEY_D))
                moveVec.x = moveSpeed * _dt;
        
            camDist -= zoomSpeed * GetMouseWheelMove();
            if (camDist < 1)
                camDist = 1;

            Vector3 camPos = { 0, 0, camDist};
            Matrix tiltMat = MatrixRotateX(tiltAngle * DEG2RAD);
            Matrix rotMat = MatrixRotateY(rotAngle * DEG2RAD); 
            Matrix mat = MatrixMultiply(tiltMat, rotMat); 

            camPos = Vector3Transform(camPos, mat);   
            moveVec = Vector3Transform(moveVec, rotMat); 
            target = Vector3Add(target, moveVec); 
            position = Vector3Add(target, camPos); 

            isDirty = false;
        }


    };

}//Util

using json = nlohmann::json;

namespace Disk {

    template<class T, template<typename> class COSTFUNCTION, class HEIGHTMAP_T, class R_T, uint32_t BUCKET_SIZE, class HEIGHTMAP_ALLOCATOR>
    TP::Config<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR> create_default() {

        using namespace TP;

        Detail::BoxList<T> postpacs;
        postpacs.shuffle_boxes_ = true;
        postpacs.list_.push_back(TP::Detail::BoxEntry<T>{ glm::vec<3, T>{(T)28., (T)17.4, (T)10.}, 100 }); 
        postpacs.list_.push_back(TP::Detail::BoxEntry<T>{ glm::vec<3, T>{(T)35.5, (T)24., (T)12.5}, 100 });
        postpacs.list_.push_back(TP::Detail::BoxEntry<T>{ glm::vec<3, T>{(T)38., (T)35., (T)16.9}, 100 });
        postpacs.list_.push_back(TP::Detail::BoxEntry<T>{ glm::vec<3, T>{(T)53.5, (T)28.5, (T)16.5}, 100 });
        postpacs.list_.push_back(TP::Detail::BoxEntry<T>{ glm::vec<3, T>{(T)39., (T)13., (T)11.}, 100 });
        postpacs.list_.push_back(TP::Detail::BoxEntry<T>{ glm::vec<3, T>{(T)55.5, (T)37., (T)6.}, 100 });

        Detail::BinInfo<T> bin;
        bin.Bounds = {80., 120.}; 
        bin.Height = 120.;

        TP::Config<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR> conf;
        conf.MultiThreading = true;
        conf.NumThreads = 4;
        conf.UseRandomSeed = true;
        conf.Seed = 12341234;
        conf.Bins.push_back(bin);
        conf.BoxType = Detail::BoxGenerationType::RANDOM;
        conf.CubeRandomBoxes = true;
        conf.LookAheadSize = 50;
        conf.EmptryTries = 0;
        conf.MaxEmptryTries = 0;
        conf.AllowOverlap = false;
        conf.MinBoxVolume = 15 * 15 * 15;
        conf.MaxBoxVolume = 30 * 30 * 30;
        conf.BoxList = std::move(postpacs);

        return conf;
    };

    template<class T, template<typename> class COSTFUNCTION, class HEIGHTMAP_T, class R_T, uint32_t BUCKET_SIZE, class HEIGHTMAP_ALLOCATOR>
    TP::Config<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR> load(const std::filesystem::path& _p) {

        using namespace TP;
        using Conf = TP::Config<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR>;

        //std::cout << _config.dump(4) << std::endl;

        if(!std::filesystem::exists("config.json")) return create_default<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR>();

        std::ifstream i(_p);
        json oo;
        i >> oo;
        
        if(!oo.contains("Config")) return create_default<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR>();

        const json& o = oo["Config"];

        Conf conf;
        conf.MultiThreading = o["MultiThreading"].template get<bool>();
        conf.NumThreads = o["NumThreads"].template get<uint32_t>();
        conf.UseRandomSeed = o["UseRandomSeed"].template get<bool>();
        conf.Seed = o["Seed"].template get<uint64_t>();
        {
            const json& bins = o["Bins"];
            for(auto it = bins.begin(); it != bins.end(); ++it){
                const json& e = *it;
                Detail::BinInfo<T> bin;
                bin.Bounds = glm::vec<2, T>{ e["Bounds"]["x"].template get<T>(), e["Bounds"]["y"].template get<T>() };
                bin.Height = e["Height"];
                conf.Bins.push_back(bin);
            }
        }
        conf.BoxType = (TP::Detail::BoxGenerationType)o["BoxType"].template get<int32_t>();
        conf.AllowOverlap = o["AllowOverlap"].template get<bool>();
        conf.EmptryTries = o["EmptryTries"].template get<uint32_t>();
        conf.MaxEmptryTries = o["MaxEmptryTries"].template get<uint32_t>();
        conf.LookAheadSize = o["LookAheadSize"].template get<uint32_t>();
        conf.AllowedPermutations = o["AllowedPermutations"].template get<uint32_t>();
        conf.CubeRandomBoxes = o["CubeRandomBoxes"].template get<bool>();
        conf.MinBoxVolume = o["MinBoxVolume"].template get<T>();
        conf.MaxBoxVolume = o["MaxBoxVolume"].template get<T>();
        {
            const json& bl = o["BoxList"];
            conf.BoxList.shuffle_boxes_ = bl["shuffle_boxes_"].template get<bool>();
            const json& list = bl["list_"];
            for(auto it = list.begin(); it != list.end(); ++it){
                const json& e = *it;
                TP::Detail::BoxEntry<T> ee;
                {
                    const json& si = e["size_"];
                    ee.size_ = glm::vec<3, T>{ si["x"].template get<T>(), si["y"].template get<T>(), si["z"].template get<T>() };
                }
                ee.perms_ = e["perms_"].template get<uint32_t>();
                ee.count_ = e["count_"].template get<uint32_t>();
                conf.BoxList.list_.push_back(ee);
            }
        }
        conf.EnforceMisses = o["EnforceMisses"].template get<bool>();
        conf.EvalBoxCount = o["EvalBoxCount"].template get<uint32_t>();
        return conf;
    }//load

    template<class T, template<typename> class CF, class H_T, class R_T, uint32_t BS, class HA>
    json save(
        const TP::Config<T, CF, H_T, R_T, BS, HA>& _conf
    ) {

        json o;
        o["MultiThreading"] = _conf.MultiThreading;
        o["NumThreads"] = _conf.NumThreads;
        o["UseRandomSeed"] = _conf.UseRandomSeed;
        o["Seed"] = _conf.Seed;
        o["Bins"] = json::array();
        for(const auto& b : _conf.Bins){
            json be = json::object();
            be["Bounds"] = { {"x", b.Bounds.x }, {"y", b.Bounds.y } };
            be["Height"] = b.Height;
            o["Bins"].push_back(be); 
        }
        o["BoxType"] = int32_t(_conf.BoxType);
        o["AllowOverlap"] = _conf.AllowOverlap;
        o["EmptryTries"] = _conf.EmptryTries;
        o["MaxEmptryTries"] = _conf.MaxEmptryTries;
        o["LookAheadSize"] = _conf.LookAheadSize;
        o["AllowedPermutations"] = _conf.AllowedPermutations;
        o["CubeRandomBoxes"] = _conf.CubeRandomBoxes;
        o["MinBoxVolume"] = _conf.MinBoxVolume;
        o["MaxBoxVolume"] = _conf.MaxBoxVolume;
        json bl = json::object();
        bl["list_"] = json::array();
        bl["shuffle_boxes_"] = _conf.BoxList.shuffle_boxes_;
        for(const TP::Detail::BoxEntry<T>& l : _conf.BoxList.list_){
            json be = json::object();
            be["size_"] = { { "x", l.size_.x }, { "y", l.size_.y }, { "z", l.size_.z } };
            be["count_"] = l.count_;
            be["perms_"] = l.perms_;
            bl["list_"].push_back(be);    
        }
        o["BoxList"] = bl;
        o["EnforceMisses"] = _conf.EnforceMisses;
        o["EvalBoxCount"] = _conf.EvalBoxCount;

        json out;
        out["Config"] = o;
        return out;
    }//save

}//Json
