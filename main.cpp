#include "iostream"
#include "../EHApplication/EHApplication.h"
#include "engine.h"
#include "fluidgroup.h"
#include <set>
#include <map>



#define MAX_FLUID_PARTICLE 1000
struct MyFrame : EH::Application::Frame
{
    FluidEngine engine;
    FluidGroup fluidgroup;
    FluidGroup staticgroup;

    std::set< std::pair< int , int > > static_map;

    void Load()
    {
        engine.Load( { -3.0f , -3.0f } , { 3.0f , 3.0f } , MAX_FLUID_PARTICLE );
        fluidgroup.particle.Load( MAX_FLUID_PARTICLE );
        engine.AddGroup( &fluidgroup );
        staticgroup.particle.Load( MAX_FLUID_PARTICLE );
        engine.AddGroup( &staticgroup );
        //grp1.gravity = { 0.0f , -5.0f };

        LoadParticleRenderProgram();
        LoadParticleMatrixProgram();

    }

    EH::GL::Program particle_render_program;
    void LoadParticleRenderProgram()
    {
        const auto m = EH::Matrix::Util::template OrthologMatrix< 4 >(
                EH::Matrix::Vector< GLfloat , 2 >( -3 ) , EH::Matrix::Vector< GLfloat , 2 >( 3 ) ,
                EH::Matrix::Vector< GLfloat , 2 >( -1 ) , EH::Matrix::Vector< GLfloat , 2 >( 1 )
                );
        m.Log();

        particle_render_program.Load();
        particle_render_program.AttachShader(
                EH::GL::LoadShaderFile( "/home/ehwan/Documents/Workspace/fluid/render_particle.vs" , GL_VERTEX_SHADER )
                );
        particle_render_program.AttachShader(
                EH::GL::LoadShaderFile( "/home/ehwan/Documents/Workspace/fluid/render_particle.fs" , GL_FRAGMENT_SHADER )
                );
        particle_render_program.BindAttribLocation( { "position" , "modelview" , "dv_white_in" } );
        particle_render_program.Link();

        particle_render_program.SetUniformMatrix4( "u_projection" , m.data() );
    }
    void RenderParticle()
    {
        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );

        particle_render_program.Use();

        {
            EH::GL::Vertex( 0 ).
                EnableVertexAttribArray().
                VertexAttribDivisor( 1 ).
                VertexAttribPointer( 2 , GL_FLOAT , 0 , fluidgroup.particle.position.get() );
            EH::GL::Vertex( 1 ).
                EnableVertexAttribArray().
                VertexAttribDivisor( 1 ).
                VertexAttribPointer( 4 , GL_FLOAT , 0 , modelview_matrix );
            EH::GL::Vertex( 2 ).
                EnableVertexAttribArray().
                VertexAttribDivisor( 1 ).
                VertexAttribPointer( 1 , GL_FLOAT , 0 , fluidgroup.particle.veldiffpressure.get() );

            glDrawArraysInstanced( GL_TRIANGLE_STRIP , 0 , 4 , fluidgroup.particle.count );

            EH::GL::Vertex( 0 ).DisableVertexAttribArray().VertexAttribDivisor();
            EH::GL::Vertex( 1 ).DisableVertexAttribArray().VertexAttribDivisor();
            EH::GL::Vertex( 2 ).DisableVertexAttribArray().VertexAttribDivisor();
        }

        glDisable( GL_BLEND );
    }

    EH::GL::Program particle_matrix_compute_program;
    EH::GL::VertexBuffer modelview_matrix;

    void LoadParticleMatrixProgram()
    {
        particle_matrix_compute_program.Load();
        particle_matrix_compute_program.AttachShader(
                EH::GL::LoadShaderFile( "/home/ehwan/Documents/Workspace/fluid/render_matrix.vs" , GL_VERTEX_SHADER )
                );
        particle_matrix_compute_program.TransformFeedbackVaryings( { "outmat" } );
        particle_matrix_compute_program.BindAttribLocation( { "velocity" , "surfacenormal" , "rho" } );
        particle_matrix_compute_program.Link();

        modelview_matrix.Load( sizeof( GLfloat )*4 * MAX_FLUID_PARTICLE , GL_DYNAMIC_DRAW , 0 );
        modelview_matrix.End();
    }
    void ParticleModelviewFeedback()
    {

        particle_matrix_compute_program.Use();

        glEnable( GL_RASTERIZER_DISCARD );
        modelview_matrix.BindTransformFeedback( 0 );

        glBeginTransformFeedback( GL_POINTS );

        EH::GL::Vertex( 0 ).EnableVertexAttribArray().VertexAttribPointer( 2 , GL_FLOAT , 0 , fluidgroup.particle.velocity.get() );
        EH::GL::Vertex( 1 ).EnableVertexAttribArray().VertexAttribPointer( 2 , GL_FLOAT , 0 , fluidgroup.particle.surfacenormal.get() );
        EH::GL::Vertex( 2 ).EnableVertexAttribArray().VertexAttribPointer( 1 , GL_FLOAT , 0 , fluidgroup.particle.rho.get() );

        glDrawArrays( GL_POINTS , 0 , fluidgroup.particle.count );

        EH::GL::Vertex( 0 ).DisableVertexAttribArray();
        EH::GL::Vertex( 1 ).DisableVertexAttribArray();
        EH::GL::Vertex( 2 ).DisableVertexAttribArray();

        glEndTransformFeedback();

        glDisable( GL_RASTERIZER_DISCARD );

        glFlush();

        EH::GL::CheckError( "ParticleModelviewFeedback" );
    }
    ~MyFrame() override
    {
        std::cout << "Destructor\n";
    }
    void EnterFrame() override
    {
        constexpr const static float HEXAGONAL_GAP = 0.1f;
        const static EH::Matrix::Matrix< GLfloat , 2 > to_hexagonal_matrix = EH::Matrix::Matrix< GLfloat , 2 >(
                EH::Matrix::Complex::Complex( 0.0f ) , EH::Matrix::Complex::Complex( (float)std::atan(1)*4.0f/3.0f )
            ) * HEXAGONAL_GAP;
        const static EH::Matrix::Matrix< GLfloat , 2 , 2 > hexagonal_coord_matrix =
            EH::Matrix::Inverse(
                to_hexagonal_matrix
            );
        if( EH::Application::Touch::IsDown( 1 ) )
        {
            const  EH::Matrix::Vector< int , 2 > hcoord = ( hexagonal_coord_matrix *  EH::Application::Touch::GetTouchPoint( 0 ) ).Convert< int >();
            if( static_map.emplace( hcoord.x , hcoord.y ).second == true )
            {
                engine.AddParticle(
                        ( to_hexagonal_matrix * hcoord ).Convert< GLfloat >() , 0 , &staticgroup );
            }
        }

        engine.Step( EH::Application::Time::GetDT() );

        if( fluidgroup.particle.count == 0 ){ return; }

        ParticleModelviewFeedback();

        /*
        modelview_matrix.Begin();
        float *mapb = (float*)modelview_matrix.MapBuffer( GL_READ_ONLY );
        EH::LOG::LOG( mapb[0] , mapb[1] , mapb[2] , mapb[3] );
        EH::LOG::LOG( mapb[4] , mapb[5] , mapb[6] , mapb[7] );
        modelview_matrix.UnmapBuffer();
        modelview_matrix.End();
        */

        glClear( GL_COLOR_BUFFER_BIT );
        EH::GL::CheckError( "glClear" );

        RenderParticle();
    }
    void TouchDown() override
    {
        const auto& p = EH::Application::Touch::GetTouchPoint();
        EH::LOG::LOG( "TouchDown ( " , EH::Application::Touch::GetID() , " ) At : " , p.x , " , " , p.y );
    }
    void TouchUp() override
    {
    }
    void TouchMove() override
    {
        if( EH::Application::Touch::IsDown( 0 ) )
        {
            engine.AddParticle( EH::Application::Touch::GetTouchPoint( 0 ) , 0.0f , &fluidgroup );
        }
    }
};

#include <memory>

void Test()
{
    using namespace EH;
    Application::Time::SetDuration(
            std::chrono::duration< float >( 1.0f/60.0f ) ,
            std::chrono::duration< float >( 1.0f/30.0f )
            );
    Application::Start();

    constexpr int w = 700;
    constexpr int h = 700;
    //Application::SetNativeResolution( { w , h } );
    Application::SetResolution( { -3.0f , -3.0f } ,  { 3.0f , 3.0f } );

    std::shared_ptr< MyFrame > frame( new MyFrame );
    Application::SetCurrentFrame( frame.get() );

    Window::InitWindow( "Test Window" , 1 , 1 , w , h );
    GL::MakeContext();
    Application::ViewportInit();
    frame->Load();
    while( EH::Window::ShouldClose() == false )
    {
        EH::Window::ProcessEvents();
        if( Application::Loop() )
        {
            Window::SwapBuffers();
        }
    }
    Window::Close();
}

int main()
{
    Test();
}
