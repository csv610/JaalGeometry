/****************************************************************************
** Meta object code from reading C++ file 'MeshVoxelizerDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshVoxelizerDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshVoxelizerDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshVoxelizerDialog_t {
    QByteArrayData data[11];
    char stringdata0[111];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshVoxelizerDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshVoxelizerDialog_t qt_meta_stringdata_JMeshVoxelizerDialog = {
    {
        QT_MOC_LITERAL(0, 0, 20), // "JMeshVoxelizerDialog"
        QT_MOC_LITERAL(1, 21, 7), // "genMesh"
        QT_MOC_LITERAL(2, 29, 0), // ""
        QT_MOC_LITERAL(3, 30, 6), // "slicer"
        QT_MOC_LITERAL(4, 37, 9), // "showEvent"
        QT_MOC_LITERAL(5, 47, 11), // "QShowEvent*"
        QT_MOC_LITERAL(6, 59, 1), // "e"
        QT_MOC_LITERAL(7, 61, 18), // "openMeshlistDialog"
        QT_MOC_LITERAL(8, 80, 6), // "saveAs"
        QT_MOC_LITERAL(9, 87, 8), // "readFile"
        QT_MOC_LITERAL(10, 96, 14) // "monotoneVoxels"

    },
    "JMeshVoxelizerDialog\0genMesh\0\0slicer\0"
    "showEvent\0QShowEvent*\0e\0openMeshlistDialog\0"
    "saveAs\0readFile\0monotoneVoxels"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshVoxelizerDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    7,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   49,    2, 0x08 /* Private */,
    3,    0,   50,    2, 0x08 /* Private */,
    4,    1,   51,    2, 0x08 /* Private */,
    7,    0,   54,    2, 0x08 /* Private */,
    8,    0,   55,    2, 0x08 /* Private */,
    9,    0,   56,    2, 0x08 /* Private */,
    10,    0,   57,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 5,    6,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMeshVoxelizerDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshVoxelizerDialog *_t = static_cast<JMeshVoxelizerDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->genMesh();
            break;
        case 1:
            _t->slicer();
            break;
        case 2:
            _t->showEvent((*reinterpret_cast< QShowEvent*(*)>(_a[1])));
            break;
        case 3:
            _t->openMeshlistDialog();
            break;
        case 4:
            _t->saveAs();
            break;
        case 5:
            _t->readFile();
            break;
        case 6:
            _t->monotoneVoxels();
            break;
        default:
            ;
        }
    }
}

const QMetaObject JMeshVoxelizerDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshVoxelizerDialog.data,
        qt_meta_data_JMeshVoxelizerDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshVoxelizerDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshVoxelizerDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshVoxelizerDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshVoxelizerDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshVoxelizerDialog"))
        return static_cast< Ui::MeshVoxelizerDialog*>(const_cast< JMeshVoxelizerDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshVoxelizerDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
