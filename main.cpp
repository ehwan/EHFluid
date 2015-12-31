#include "iostream"
#include "../EHApplication/EHApplication.h"
#include "engine.h"
#include "fluidgroup.h"


struct MyFrame : EH::Application::Frame
{
    FluidEngine engine;
    FluidGroup grp1;

    EH::GL::Program program;

    void Load()
    {
        engine.Load( { -3.0f , -3.0f } , { 3.0f , 3.0f } , 1000 );
        grp1.particle.Load( 1000 );
        engine.AddGroup( &grp1 );
        //grp1.gravity = { 0.0f , -5.0f };

        const char * vs=
            "#version 130\n"
            "in vec2 vert;\n"
            "void main()"
            "{"
            "gl_Position = vec4( vert/3.0 , 0.0 , 1.0 );"
            "}";
        const char *fs =
            "#version 130\n"
            "void main()"
            "{"
            "gl_FragColor=vec4(1.0);"
            "}";

        program.Load();
        program.AttachShader( EH::GL::LoadShader( vs , GL_VERTEX_SHADER ) );
        program.AttachShader( EH::GL::LoadShader( fs , GL_FRAGMENT_SHADER ) );
        glBindAttribLocation( program.program , 0 , "vert" );
        program.Link();
        ::glPointSize( 4.0f );
    }
    ~MyFrame() override
    {
        std::cout << "Destructor\n";
    }
    void EnterFrame() override
    {
        engine.Step( EH::Application::Time::GetDT() );

        glClear( GL_COLOR_BUFFER_BIT );
        program.Use();
        ::glEnableVertexAttribArray( 0 );
        ::glVertexAttribPointer( 0 , 2 , GL_FLOAT , GL_FALSE , 0 , grp1.particle.position.get() );
        ::glDrawArrays( GL_POINTS , 0 , grp1.particle.count );
        ::glDisableVertexAttribArray( 0 );
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
            engine.AddParticle( EH::Application::Touch::GetTouchPoint( 0 ) , 0.0f , &grp1 );
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
